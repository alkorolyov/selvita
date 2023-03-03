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


def colorize(substring, text):
    matches = find_near_matches(substring, text, max_l_dist=2)
    if matches:
        return text.replace(matches[0].matched, f'{Fore.BLACK}{Back.LIGHTYELLOW_EX}{substring}{Style.RESET_ALL}')
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


def draw_reaction(data: Union[pd.Series, pd.DataFrame], highlight_text: str = None, render_format='png'):
    if isinstance(data, pd.DataFrame):
        data = data.sample().squeeze()

    rxn = indigo.loadReaction(data['rxn_smiles'])
    indigo.setOption("render-output-format", render_format)
    if render_format == 'png':
        display(Image(renderer.renderToBuffer(rxn)))
    elif render_format == 'svg':
        display(SVG(renderer.renderToBuffer(rxn)))
    else:
        print(f"{render_format}: unknown render format, 'svg' or 'png' ")

    print("Patent:      ", data['patent'])
    print("Reaction_id: ", data.name)

    if highlight_text:
        print(colorize(highlight_text, data['notes']))
    else:
        print(data['notes'])


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
