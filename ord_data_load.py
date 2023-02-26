import os
import re
import json
import gzip
import multiprocessing as mp

import pandas as pd
import numpy as np
from google import protobuf
from ord_schema import message_helpers
from ord_schema.proto import dataset_pb2

from copy import deepcopy

from rdkit.Chem import MolFromSmiles as smiles2mol
from rdkit.Chem.AllChem import ReactionFromSmarts
from rdkit.Chem import MolFromSmiles
from rdkit.Chem.Draw import MolsToGridImage
from rdkit.Chem import Draw
from rdkit.Chem.Draw import rdMolDraw2D

from indigo import Indigo
from indigo.renderer import IndigoRenderer

from IPython.display import SVG
from IPython.display import display_svg
from IPython.display import display_png
from IPython.display import Image
from IPython.display import display

from typing import Union

indigo = Indigo()
renderer = IndigoRenderer(indigo)

indigo.setOption("render-output-format", "svg")
# indigo.setOption("render-superatom-mode", "collapse")
indigo.setOption("render-coloring", True)
# indigo.setOption("render-base-color", "1, 1, 1")
indigo.setOption("render-relative-thickness", 1.5)


ORD_REPO_PATH = './ord-data'
ORD_PATH = './ORD'


class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKCYAN = '\033[96m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'


def colorize(text, substring, color=bcolors.WARNING):
    return text.replace(substring, f'{color}{substring}{bcolors.ENDC}')


def df_na_vals(df, return_empty=True):
    columns = df.columns
    N = max(len(c) for c in columns) + 5
    empty = []
    for col in columns:
        na_vals = df[col].isna() | ~df[col].apply(bool)
        print(col.ljust(N), '->', ' '*(N//3), f'Missing values: {na_vals.sum()} ({na_vals.mean():.2%})')

        if na_vals.mean() > .99:
            empty.append(col)

    if return_empty:
        return empty


def fahreneheit_to_celsius(t):
    return (t - 32) * (5/9)


def extract_temperature_from_notes(description, room_temperature=25):
    if not description:
        return None
    match = re.search('\s(?P<temperature>-{0,1}\d{1,3})[\s°o\*]*(c|C|deg|degrees)\W', description)
    if match:
        return int(match.groupdict()['temperature'])
    elif re.search('\W(rt|RT|room temp|ambient temp)', description):
        return room_temperature
    else:
        return None


def extract_reaction_conditions(dic, roles_map={'REACTANT': 'reactants', 'SOLVENT': 'solvents', 'CATALYST': 'catalysts'}):
    # Reaction identifier
    final_dic = {'id': dic['reactionId']}

    # Reaction SMILES
    for i in dic['identifiers']:
        if i['type'] == 'REACTION_CXSMILES':
            final_dic['reaction_smile'] = i['value'].split()[0] # remove the
            break
    else:
        final_dic['reaction_smile'] = None

    # Compounds
    final_dic['solvents'] = []
    final_dic['reactants'] = []
    final_dic['catalysts'] = []
    final_dic['reagents'] = []
    for number, input_dics in dic['inputs'].items():
        for i_d in input_dics['components']:
            role = roles_map.get(i_d['reactionRole'], 'reagents')
            compound_name = None
            compound_smiles = None
            for i in i_d['identifiers']:
                if i['type'] == 'NAME':
                    compound_name = i['value'].replace('′',"'")
                elif i['type'] == 'SMILES':
                    compound_smiles = i['value']
            final_dic[role].append((compound_name, compound_smiles))

    # Temperature
    temperature = dic['conditions'].get('temperature', [])
    if 'control' in temperature:
        temperature = temperature['control']['type']
    elif 'setpoint' in temperature:
        t = temperature['setpoint']['value']
        if temperature['setpoint']['units'] == 'KELVIN':
            temperature = t - 273.15
        elif temperature['setpoint']['units'] == 'FAHRENHEIT':
            temperature = fahreneheit_to_celsius(t)
        else:
            temperature = t
    else:
        temperature = extract_temperature_from_notes(dic['notes'].get('procedureDetails', ''), 'AMBIENT')

    final_dic['temperature'] = temperature


    # Time
    outcome = dic['outcomes'][0]
    time = outcome.get('reactionTime')
    if time:
        if time['units'] == 'MINUTES':
            time = time['value'] / 60
        else:
            time = time['value']

    final_dic['time'] = time

    # Product and Yield
    for product in outcome['products']:
        for p_id in product.get('identifiers', []):
            if p_id['type'] == 'SMILES':
                product_smiles = p_id['value']
                break
        else:
            product_smiles = None

        measurements = product.get('measurements', [])
        for measurement in measurements:
            if measurement['type'] == 'YIELD':
                reaction_yield = measurement['percentage']['value']
                break
        else:
            reaction_yield = None

    final_dic['product'] = product_smiles
    final_dic['yield'] = reaction_yield

    # Notes
    final_dic['notes'] = dic.get('notes', {}).get('procedureDetails')

    # Provenance
    final_dic['patent'] = dic.get('provenance', {}).get('patent')

    return final_dic


def parse_pb_file(pb: str, ord_parsed_path: str):
    # make dirs
    os.makedirs(f'{ord_parsed_path}/originals/', exist_ok=True)
    os.makedirs(f'{ord_parsed_path}/parsed/', exist_ok=True)
    dataset_name = os.path.split(pb)[-1].strip('.pb.gz')
    reaction_conditions = []
    try:
        data = message_helpers.load_message(pb, dataset_pb2.Dataset)
        data_dic = message_helpers.json_format.MessageToDict(data)
        with open(f'{ord_parsed_path}/originals/{dataset_name}.json', 'w') as file:
            json.dump(data_dic, file)
        for rxn_dic in data_dic['reactions']:
            reaction_conditions.append(extract_reaction_conditions(rxn_dic))
        # with mp.Pool(n_cores) as p:
        #     reaction_conditions = p.map(extract_reaction_conditions, data_dic['reactions'])
    except Exception as e:
        print(f'{dataset_name} extraction failed with error: {e}')

    if reaction_conditions:
        with open(f'{ord_parsed_path}/parsed/{dataset_name}.json', 'w') as file:
            json.dump(reaction_conditions, file)


def load_dataset(filename: str) -> dataset_pb2.Dataset:
    try:
        with gzip.open(filename, "rb") as f:
            return dataset_pb2.Dataset.FromString(f.read())
    except protobuf.message.DecodeError as error:
        raise ValueError(f"error parsing {filename}: {error}") from error


def filter_uspto_filenames(filename: str) -> Union[str, None]:
    try:
        dataset = load_dataset(filename)
        if "uspto" in dataset.name:
            return filename
    except protobuf.message.DecodeError as error:
        raise ValueError(f"error parsing {filename}: {error}") from error


def is_reaction_of_type(reaction_to_test,
                        # reaction_type_pattern=ReactionFromSmarts("[#8]-[#5](-[#8])-[#6:1].[#17,#35,#53]-[#6:2]>>[#6:1]-[#6:2]"),
                        reaction_type_pattern=None,
                        *reactants_patterns,
                        max_products=10
                        ):
    if reaction_to_test is None:
        return False
    if isinstance(reaction_to_test, str):
        reactants = tuple(filter(None, (smiles2mol(smiles) for smiles in reaction_to_test.split('>')[0].split('.'))))
        actual_products = tuple(filter(None, (smiles2mol(smiles) for smiles in reaction_to_test.split('>')[-1].split('.'))))
    else:
        reactants = reaction_to_test.GetReactants()
        actual_products = reaction_to_test.GetProducts()

    if isinstance(reaction_type_pattern, str):
        # use input string to make an RdKit reaction and molecules
        reactants_patterns = tuple(smiles2mol(smiles) for smiles in reaction_type_pattern.split('>')[0].split('.'))
        reaction_type_pattern = ReactionFromSmarts(reaction_type_pattern)

    elif len(reactants_patterns) == 0:
        # target reactants to be extracted from the target reaction
        reactants_patterns = reaction_type_pattern.GetReactants()

    elif isinstance(reactants_patterns[0], str):
        # target reactants to be cast as RdKit Molecules
        reactants_patterns = tuple(smiles2mol(smiles) for smiles in reactants_patterns)

    else:
        # everything is already an RdKit object
        pass

    reactants_to_test = []
    for r_pattern in reactants_patterns:
        for r, reactant in enumerate(reactants):
            if reactant is None:
                continue
            reactant.UpdatePropertyCache()
            if reactant.HasSubstructMatch(r_pattern):
                # the reactant matches with the pattern
                reactants_to_test.append(deepcopy(reactant))
                break # stop searching
        else:
            # no reactant in the input reaction matched with the pattern
            return False

    for p in actual_products:
        p.UpdatePropertyCache()

    estimated_products = []
    for product, *_ in reaction_type_pattern.RunReactants(reactants_to_test, maxProducts=max_products):
        try:
            product.UpdatePropertyCache()
            estimated_products.append(product)
        except:
            pass

    if len(estimated_products) == 0:
        return False

    return any(any(pred.HasSubstructMatch(obs) for obs in actual_products)  # any of the actual products matches the predicted product
               for pred in estimated_products)  # for any of the predicted products


def clear_atom_mapping(rxn):
    for m in rxn.GetReactants():
        for atom in m.GetAtoms():
            if atom.GetAtomMapNum() != 0:
                atom.ClearProp('molAtomMapNumber')
    for m in rxn.GetProducts():
        for atom in m.GetAtoms():
            if atom.GetAtomMapNum() != 0:
                atom.ClearProp('molAtomMapNumber')


def draw_reaction_solvents(df):
    rxn_index = np.random.randint(0, df.shape[0])
    rxn = df.loc[df.index[rxn_index]]
    rxn_mol = ReactionFromSmarts(rxn['reaction_smile'], useSmiles=True)
    clear_atom_mapping(rxn_mol)
    do = rdMolDraw2D.MolDrawOptions()
    do.includeAtomTags = False
    do.includeMetadata = False
    print('Solvent:', rxn['solvents'])
    display(Draw.ReactionToImage(rxn_mol,
                                 subImgSize=(300, 300),
                                 drawOptions=do))
    notes = rxn['notes']
    for (name, _) in rxn['solvents']:
        notes = colorize(notes, name)
    print(notes)


def draw_reaction(data: Union[pd.Series, pd.DataFrame], highlight_text: str = None, render_format='png'):
    if isinstance(data, pd.DataFrame):
        data = data.sample()

    rxn = indigo.loadReaction(data['reaction_smile'].item())
    indigo.setOption("render-output-format", render_format)
    if render_format == 'png':
        display(Image(renderer.renderToBuffer(rxn)))
    elif render_format == 'svg':
        display(SVG(renderer.renderToBuffer(rxn)))
    else:
        print(f"{render_format}: unknown render format, 'svg' or 'png' ")

    print("Patent:      ", data['patent'].item())
    print("Reaction_id: ", data.index.item())

    if highlight_text:
        print(colorize(data['notes'].item(), highlight_text))
    else:
        print(data['notes'].item())

def draw_reaction_smi(rxn_smiles: str, render_format='png'):
    rxn = indigo.loadReaction(rxn_smiles)
    indigo.setOption("render-output-format", render_format)
    if render_format == 'png':
        display(Image(renderer.renderToBuffer(rxn)))
    elif render_format == 'svg':
        display(SVG(renderer.renderToBuffer(rxn)))
    else:
        print(f"{render_format}: unknown render format, 'svg' or 'png' ")

def draw_mol(mol_smi: str, render_format='png'):
    mol = indigo.loadMolecule(mol_smi)
    indigo.setOption("render-output-format", render_format)
    if render_format == 'png':
        display(Image(renderer.renderToBuffer(mol)))
    elif render_format == 'svg':
        display(SVG(renderer.renderToBuffer(mol)))
    else:
        print(f"{render_format}: unknown render format, 'svg' or 'png' ")