# -*- coding: utf-8 -*-

#  This script was made by CCD
#  struc_type_v1.0_20221020
#  Condensed Matter Physics Group

#import module
from ase import Atom,Atoms
from ase.constraints import FixAtoms
from ase.io.vasp import write_vasp,read_vasp
from ase.build import sort, add_vacuum
from random import shuffle
from pymatgen.core import Lattice, Structure, Molecule
import numpy as np
import math
import subprocess
import re

def list_del(list_to_be_removed, list_contain_remove_item):
	for i in list_contain_remove_item:								#remove all structure saved to dict 
		list_to_be_removed.remove(i)
	return list_to_be_removed
	
def top_bot_sym(struc1):
#read CONTCAR
	CONTCAR = read_vasp(struc1)
	cell = CONTCAR.get_cell()
	#print("CONTCAR:", CONTCAR)
	#print("cell", cell)

#remove copper surface
	rm_list = []
	for atom in CONTCAR:
		if atom.position[2] < 0.05 :
			rm_list.append(atom.index)
	#print("Cu surface index:", rm_list)
	del CONTCAR[rm_list]
	#print("CONTCAR rm Cu surface:", CONTCAR)

#identify top and bottom layers
	top_layer = []
	bottom_layer = []
	sort_list = np.sort(CONTCAR.get_positions()[:,2])
	#print("sorted z-axis pos:", sort_list)

	for atom in CONTCAR:
		if abs(atom.position[2] - sort_list[0]) < 0.5:
			bottom_layer.append(atom.index)
	
	for atom in CONTCAR:
		if abs(atom.position[2] - sort_list[-1]) < 0.5:
			top_layer.append(atom.index)

	#print("bot_layer index:", bottom_layer)
	#print("top_layer index:", top_layer)

#confirm that the number of atoms in each layer is the same
	if len(bottom_layer) != len(top_layer):
		#print("atom numbers aren't equal!", "bottom layer has", len(bottom_layer), "atoms not equal to top layer has", len(top_layer), "atoms")
		return "not sym"
	else:
		#print("the number of atoms in each layer is the same!")

#confirm position
		for index1 in bottom_layer:
			#print("bot index", index1)
			pos1 = CONTCAR.get_scaled_positions()[index1]
			match_time = 0
			#print("bot pos:", pos1)
			for index2 in top_layer:
				#print("top index", index2)
				pos2 = CONTCAR.get_scaled_positions()[index2]
				#print("top pos:", pos2)
				vec = abs(pos1 - pos2)
				#print("vec diff:", vec)
				for i in range(3):	
					if vec[i] > 0.5:
						vec[i] -=1
				#print("vec diff wrap:", vec)
				dot = np.dot(vec, cell)							#transform to Cartesian
				distance = np.linalg.norm(dot[0:-1])					#calculate 2d(x,y) distance
				#print("2d distance:", distance)

#check distance
				if distance < 0.2:
#check atom's symbol are equal
					if CONTCAR[index1].symbol == CONTCAR[index2].symbol:
						match_time += 1
						#print("match")
						top_layer.remove(index2)				#remove matched atom
						#print("top_layer index after remove:", top_layer)
						break
					#else:
						#print("atom symbols aren't equal!", CONTCAR[index1].symbol, "not equal to", CONTCAR[index2].symbol)

#check atom has been matched
			if match_time == 0:
				#print("index", index1, "not matched")
				return "not sym"
				break
		return "sym"

def struc_match_in_list(folder_list, path = "/OPTCELL/CONTCAR", file_name = 'structure_type'):
	
#creat dict to save different structure list
	struc_dict = {}

	while True:
		print(len(folder_list), "struc in list")
		if len(folder_list) == 0:							#break while loop when list is empty
			print("match completed\n")
			break
		comparison_struc = Structure.from_file(folder_list[0] + path)			#choose index1 in folder as comparison structure
		#print("comparison_structure", comparison_struc)
		same_list = []									#creat list to save the same structure as the comparison structure
		for i in folder_list:
			struc = Structure.from_file(i + path)
			if comparison_struc.matches(struc) == True:				#matches
				same_list.append(i)						#add the same structure to list
		print(len(same_list), "same struc:\n", same_list, "\n")
		struc_dict[folder_list[0]] = same_list						#add list to dict
		
		list_del(folder_list, same_list)
		#print("after remove", folder_list)
	#print("dict", struc_dict)

#clear file 
	with open(file_name, 'w') as f:
		print("clear file:", file_name)

#write file
	print("\n" + str(len(struc_dict)), "struc type:")
	for k in struc_dict:
		print(re.sub("\'|\,|\[|\]", "", str(struc_dict[k])))
		#print(re.sub("\'", "", struc_dict[k]))
		with open(file_name, 'a') as f:
			print(re.sub("\'|\,|\[|\]", "", str(struc_dict[k])), file = f)	
