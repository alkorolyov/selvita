import os
import re
import json
import gzip

import numpy as np
import pandas as pd

from google import protobuf
from ord_schema import message_helpers
from ord_schema.proto import dataset_pb2
from ord_schema.proto import reaction_pb2

from copy import deepcopy

from rdkit.Chem import MolFromSmiles as smiles2mol
from rdkit.Chem.AllChem import ReactionFromSmarts
from rdkit.Chem.rdChemReactions import ChemicalReaction, HasReactionSubstructMatch


from typing import Union

ORD_REPO_PATH = './ord-data'
ORD_PATH = './ORD'

TEMP_CONTROL_MAP = {
    2: 25.0,        # AMBIENT
    6: 0.0,         # ICE_BATH
    9: -78.0,       # DRY_ICE_BATH
    11: -120.0      # LIQUID_NITROGEN
}


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


def fahrenheit_to_celsius(t):
    return (t - 32) * (5/9)


"""

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
            temperature = fahrenheit_to_celsius(t)
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

"""


""" NEW FORMAT OF PARSED PB FILES """



def load_dataset(filename: str) -> dataset_pb2.Dataset:
    try:
        with gzip.open(filename, "rb") as f:
            return dataset_pb2.Dataset.FromString(f.read())
    except protobuf.message.DecodeError as error:
        raise ValueError(f"error parsing {filename}: {error}") from error


def parse_compound(arr: np.array,
                   idx: int,
                   rxn: reaction_pb2.Reaction,
                   cmpd: Union[reaction_pb2.Compound, reaction_pb2.ProductCompound]
                   ):
    """
    In USPTO scheme there are two or single name identifiers:
        [systematic]
        [trivial, systematic]

        arr[:, 0] <- trivial
        arr[:, 1] <- systematic
    """
    names = []
    for i in cmpd.identifiers:
        if i.type == reaction_pb2.CompoundIdentifier.NAME:
            names.append(i.value)
        if i.type == reaction_pb2.CompoundIdentifier.SMILES:
            arr[idx, 2] = i.value

    if len(names) == 1:
        # [systematic]
        arr[idx, 1] = names[0]
    else:
        # [trivial, systematic]
        arr[idx, 0] = names[0]
        arr[idx, 1] = names[1]

    arr[idx, 3] = cmpd.reaction_role
    arr[idx, 4] = rxn.reaction_id

def parse_product(arr: np.array,
                  idx: int,
                  rxn: reaction_pb2.Reaction,
                  cmpd: Union[reaction_pb2.Compound, reaction_pb2.ProductCompound]
                  ):

    parse_compound(arr, idx, rxn, cmpd)
    # products sometimes have UNDEFINED rxn_role
    if arr[idx, 3] == 0: # UNDEFINED
        arr[idx, 3] = 8 # PRODUCT


def parse_dataset(dataset: dataset_pb2.Dataset) -> np.ndarray:
    """
        Parses USPTO dataset to extract compound data. Creates specific
        numpy array of compounds with two indexes [0 .. total_cmpd_numer, field_idx]
    - first idx:
        Compound index
    - second idx:
        Corresponding field index
        ['second name', 'name', 'smiles', 'role', 'rxn_id']
        arr[i, 0] - (str/optional), trivial name or compound label in the patent
        arr[i, 1] - (str), systematic name
        arr[i, 2] - (str), smiles
        arr[i, 3] - (int), reaction role enum from reaction_pb2.ReactionRole.ReactionRoleType,
                           e.g. "REACTANT" - 1, "SOLVENT" - 3, "CATALYST" - 4, "PRODUCT" - 8
        arr[i, 4] - (str), reaction_id, e.g. "ord-43d5b7a6265d46a0ab8a7e2b2db5ad33"

    :param dataset:
    :return:
    """
    N = len(dataset.reactions)
    arr = np.empty((N*10, 5), dtype=object)
    idx = 0

    for rxn in dataset.reactions:
        compounds = []
        products = []

        for key in rxn.inputs:
            compounds.extend(rxn.inputs[key].components)

        products.extend(rxn.outcomes[0].products)

        for cmpd in compounds:
            parse_compound(arr, idx, rxn, cmpd)
            idx += 1
        for cmpd in products:
            parse_product(arr, idx, rxn, cmpd)
            idx += 1
    arr = arr[:idx]
    return arr


def pb2_to_numpy_cmpd(filename: str) -> np.ndarray:
    """
    Helper function for multiprocessing
    """
    dataset = load_dataset(filename)
    return parse_dataset(dataset)


def parse_dataset_rxn(dataset: dataset_pb2.Dataset) -> np.ndarray:
    """
    Parses USPTO dataset to extract reaction data. Creates specific
    numpy array of compounds with two indexes [0 .. N, field_idx]

    columns = ['rxn_id', 'rxn_smiles', 'time_unit', 'time_val', 'temp_unit', 'temp_val', 'temp_control', 'yield', 'patent', 'notes']
    field_idx:   0            1             2             3            4             5              6           7        8         9

    :param dataset: ord USPTO dataset
    :return: numpy array with shape (N, 7)
    """
    N = len(dataset.reactions)
    arr = np.empty((N, 10), dtype=object)

    for idx, rxn in enumerate(dataset.reactions):
        # fields always present in USPTO
        arr[idx, 0] = rxn.reaction_id
        arr[idx, 1] = rxn.identifiers[0].value  # rxn_smiles

        # time (hours)
        if rxn.outcomes[0].HasField('reaction_time'):
            time = rxn.outcomes[0].reaction_time
            arr[idx, 2] = time.units
            arr[idx, 3] = time.value

        # temperature (°C)
        temp = rxn.conditions.temperature
        if temp.HasField('setpoint'):
            arr[idx, 4] = temp.setpoint.units
            arr[idx, 5] = temp.setpoint.value

        # ambient temp control
        if temp.HasField('control'):
            arr[idx, 6] = temp.control.type


        # yield (keep PERCENTYIELD if more than two)
        yields = {}
        for p in rxn.outcomes[0].products:
            for m in p.measurements:
                if m.type == 3:  # YIELD
                    yields[m.details] = m.percentage.value
        if yields:
            if len(yields) > 1:
                y = yields.get("PERCENTYIELD", None)
            else:
                y = yields.get("CALCULATEDPERCENTYIELD", None)
            arr[idx, 7] = y

        arr[idx, 8] = rxn.provenance.patent
        arr[idx, 9] = rxn.notes.procedure_details


    arr = arr[:idx]
    return arr


def pb2_to_numpy_rxn(filename: str) -> np.ndarray:
    """
    Helper function for multiprocessing
    """
    dataset = load_dataset(filename)
    return parse_dataset_rxn(dataset)


def pb2_test_parse(filename: str) -> np.ndarray:
    """
    Explore values and units
    :param filename:
    :return:
    """
    dataset = load_dataset(filename)
    N = len(dataset.reactions)
    arr = np.empty((N, 8), dtype=object)

    for idx, rxn in enumerate(dataset.reactions):
        if rxn.outcomes[0].HasField('reaction_time'):
            time = rxn.outcomes[0].reaction_time
            t = None
            if time.units == 1:     # HOUR
                t = time.value
            elif time.units == 2:   # MINUTE
                t = time.value / 60
            elif time.units == 3:   # SECOND
                t = time.value / 3600
            elif time.units == 4:   # DAY
                t = time.value * 24
            arr[idx, 1] = time.units
            arr[idx, 2] = t

        temp = rxn.conditions.temperature
        if temp.HasField('setpoint'):
            arr[idx, 3] = temp.setpoint.units
            arr[idx, 4] = temp.setpoint.value

        if temp.HasField('control'):
            arr[idx, 5] = temp.control.type

        arr[idx,0] = rxn.reaction_id
        arr[idx, 6] = rxn.notes.procedure_details
        arr[idx, 7] = rxn.identifiers[0].value

    arr = arr[:idx]

    return arr

def filter_uspto_filenames(filename: str) -> Union[str, None]:
    """
    Helper function for multiprocessing
    """
    try:
        dataset = load_dataset(filename)
        if "uspto" in dataset.name:
            return filename
    except protobuf.message.DecodeError as error:
        raise ValueError(f"error parsing {filename}: {error}") from error


""" REACTION SUBSTRUCTURE SEARCH """

def is_reaction_of_type(reaction_to_test: str,
                        reaction_type_pattern: ChemicalReaction=None,
                        max_products=10
                        ):
    if reaction_to_test is None:
        return False

    reactants = tuple(filter(None, (smiles2mol(smiles) for smiles in reaction_to_test.split(' ')[0].split('>')[0].split('.'))))
    actual_products = tuple(filter(None, (smiles2mol(smiles) for smiles in reaction_to_test.split(' ')[0].split('>')[-1].split('.'))))

    reactants_patterns = reaction_type_pattern.GetReactants()

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

    if not estimated_products:
        return False

    return any(any(pred.HasSubstructMatch(obs) for obs in actual_products)  # any of the actual products matches the predicted product
               for pred in estimated_products)  # for any of the predicted products