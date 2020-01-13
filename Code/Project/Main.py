from Structure import *
from JsonRead import *
from Material import *

structure_js = readFile('testStructure.json')
# structure_js = readFile('structure00.json')
material_models_js = readFile('material_models.json')
load_material_models(material_models_js)
structure = Structure(structure_js)
# structure.analyze_structure()

