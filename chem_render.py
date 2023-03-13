from typing import Union
from colorama import Fore, Back, Style
from fuzzysearch import find_near_matches

import numpy as np
import pandas as pd
from IPython.core.display import display, Image, SVG
from indigo import Indigo
from indigo.renderer import IndigoRenderer
from rdkit.Chem import Draw
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit.Chem.rdChemReactions import ReactionFromSmarts

indigo = Indigo()
renderer = IndigoRenderer(indigo)

indigo.setOption("render-output-format", "png")
# indigo.setOption("render-superatom-mode", "collapse")
indigo.setOption("render-coloring", True)
# indigo.setOption("render-base-color", "1, 1, 1")
indigo.setOption("render-relative-thickness", 1.5)
indigo.setOption("render-highlight-thickness-enabled", True)
indigo.setOption("render-highlighted-labels-visible", True)
# indigo.setOption("render-highlight-color", "1, 0.4, 0")


def colorize(substring, text):
    matches = find_near_matches(substring, text, max_l_dist=1)
    if matches:
        for m in matches:
            text = text.replace(m.matched, f'{Fore.BLACK}{Back.LIGHTYELLOW_EX}{m.matched}{Style.RESET_ALL}')
    return text


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
        notes = colorize(name, notes)
    print(notes)


def draw_reaction(data: Union[pd.Series, pd.DataFrame],
                  highlight_text: str = None,
                  highlight_smi: str = None,
                  render_format='png',
                  auto_map: bool = False) -> Union[str, None]:
    """
    Sample and draw random reaction from dataset
    :param auto_map: render atom to atom mapping by Indigo
    :param data: input dataframe
    :param highlight_text: highlighted text in reaction notes, with up to
                           1 difference in levenshtein distance.
    :param render_format: 'png' or 'svg'
    :return: sampled reaction id
    """
    if isinstance(data, pd.DataFrame):
        size = len(data)
        data = data.sample().squeeze()

    print("Set size:        ", size)
    print("Patent:          ", data.get('patent'))
    print("Reaction_id:     ", data.name)

    try:
        rxn = indigo.loadReaction(data['rxn_smiles'])
        if auto_map:
            rxn.automap()
    except Exception as e:
        print(f"Parsing error: {e}\nSmiles: {data['rxn_smiles']}")
        return None

    if highlight_smi:
        indigo.setOption("render-coloring", False)

        query = indigo.loadReactionSmarts(highlight_smi)
        match = indigo.substructureMatcher(rxn).match(query)
        if match:
            draw_indigo_obj(match.highlightedTarget())
        else:
            print(f"Pattern {highlight_smi} not found")
            draw_indigo_obj(rxn)

        indigo.setOption("render-coloring", True)

    else:
        draw_indigo_obj(rxn)

    # print("Reaction SMARTS: ", data.get('rxn_smiles'))

    if highlight_text:
        print(colorize(highlight_text, data.get('notes')))
    else:
        print(data.get('notes'))
    return data.name


def draw_indigo_obj(obj, render_format='png'):
    indigo.setOption("render-output-format", render_format)
    if render_format == 'png':
        display(Image(renderer.renderToBuffer(obj)))
    elif render_format == 'svg':
        display(SVG(renderer.renderToBuffer(obj)))
    else:
        print(f"{render_format}: unknown render format, 'svg' or 'png' ")


def draw_reaction_smi(rxn_smiles: str, render_format='png'):
    rxn = indigo.loadReaction(rxn_smiles)
    draw_indigo_obj(rxn)


def draw_mol(mol_smi: str, render_format='png'):
    mol = indigo.loadMolecule(mol_smi)
    draw_indigo_obj(mol)
