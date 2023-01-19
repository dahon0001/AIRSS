# -*- coding: utf-8 -*-

#  This script was made by CCD
#  AIRSS_v4.0_20220816
#  Condensed Matter Physics Group

#載入函數
from ase import Atom,Atoms
from ase.constraints import FixAtoms
from ase.io.vasp import write_vasp,read_vasp
from ase.build import sort, add_vacuum
from random import shuffle
import numpy as np
import math
import subprocess

#繼承Atoms
class Atoms_random(Atoms):

	def __init__(self, numbers_of_structure, append_elements, numbers_of_elements, cell):
		self.numbers_of_structure = numbers_of_structure
		self.append_elements = append_elements
		self.numbers_of_elements = numbers_of_elements
		self.cell = cell

	#兩direct位置距離
	def get_distance(positions1, positions2, cell):
		vec = abs(positions1 - positions2)						#兩位置向量差
		for i in range(3):
			if vec[i] > 0.5:
				vec[i] -= 1
		dot = np.dot(vec, cell)							#轉成Cartesian
		distance = np.linalg.norm(dot)						#計算距離
		return distance

	#確認list中所有位置距離
	def check_distance(list1, cell):
		for pos1 in range(len(list1)):
			for pos2 in range(len(list1)):
				if pos1 >= pos2:
					continue					#不用重複判斷(判斷過1對2，就不用判斷2對1)
				distance = Atoms_random.get_distance(list1[pos1], list1[pos2], cell)
				print(pos1, pos2, "distance", distance)

	#隨機結構
	def random(numbers_of_structure, append_elements, numbers_of_elements, ini_cell = [5, 5, 5], min_distance = 2, max_distance = 5, add_vac = 5, final_vac = 1.5):

		#建幾個結構
		for number in range(numbers_of_structure):
			print(number + 1, "start")

			#將元素與個數寫進list
			append_list = []
			for append in range(len(append_elements)):
				append_list += [append_elements[append]] * numbers_of_elements[append]
			shuffle(append_list)

			#建構結構
			rs = Atoms(append_list[0] + append_list[1], scaled_positions = [(0.0, 0.0, 0.0), (0.5, 0.5, 0.0)], cell = [3.6153800488, 3.6153800488, 5])

			#加入原子與位置
			for elements in append_list[2:]:
				condition = 0						#重置while True
				times = 0

				#生成亂數並判斷距離
				while condition == 0:					#無限迴圈
					random_position = np.random.rand(3)			#生成[1*3]亂數(位置)
					max_dis_times = 0				#小於最大距離次數
					times +=1
					#print(elements, times, "times random")
							
					#判斷距離
					for r_pos in range(len(rs.get_scaled_positions())):
						distance = Atoms_random.get_distance(random_position, rs.get_scaled_positions()[r_pos], rs.get_cell())
						#print(r_pos, "r_pos:", rs.get_scaled_positions()[r_pos], "random:", random, "distance:", distance)
							
						#距離太近
						if distance <= min_distance:
							#print(distance, "< min, retry")
							break				#打斷for迴圈，回到while迴圈重新生成亂數

						#距離太遠
						if distance <= max_distance:		#判斷小於最大距離
							max_dis_times += 1		#次數+1
							#print(distance, "< max", max_dis_times, "times")
					
						#確認與所有原子距離合適
						if r_pos == len(rs.get_scaled_positions()) - 1:			#確認所有位置距離都判斷過
							#print("all positions are checked")

							if max_dis_times > 0:									#確認最近原子距離在最大距離內
								rs.append(elements)
								rs[-1].scaled_position = random_position
								rs.center(vacuum = 0, axis = 2)									#將cluster移到中心並加上真空層
								add_vacuum(rs, 5)							#將cluster移到中心並加上真空層
								#rs = sort(rs)									#sort symbol
								condition += 1									#使while False，打斷迴圈
								#write_vasp(str(rs.symbols) + "_" + "{:02d}".format(number + 1) + ".vasp", rs, direct=True, sort=True, vasp5=True)
								#print("ok, < max", max_dis_times, "times")

							else:				#最近原子不在最大距離內
								#print("retry, < max", max_dis_times, "times")
								break			#打斷for迴圈，回到while迴圈重新生成亂數
				
			#確認距離
			#Atoms_random.check_distance(rs.get_scaled_positions(), rs.get_cell())

			#調整真空層
			rs.center(vacuum = 0, axis = 2)
			add_vacuum(rs, 3)
			top_layer = Atoms(rs.symbols[0] + rs.symbols[1], scaled_positions = [(0.0, 0.0, 1.0), (0.5, 0.5, 1.0)], cell= rs.get_cell())
			rs += top_layer
			rs.center(vacuum = 3, axis=2)
			Cu = Atoms("Cu2", scaled_positions = [(0.5, 0.0, 0.0), (0.0, 0.5, 0.0)], cell = rs.get_cell())
			rs += Cu
			rs.set_constraint(FixAtoms(indices=[atom.index for atom in rs if atom.position[2] < 0.05]))
			rs = sort(rs)


			#輸出vasp檔
			subprocess.run("mkdir " + str(rs.symbols) + "_" + "{:02d}".format(number + 1), shell = True)	#創建資料夾
			write_vasp(str(rs.symbols) + "_" + "{:02d}".format(number + 1) + "/POSCAR_0", rs, direct=True, sort=True, vasp5=True)
			print(number + 1, "complete")

Atoms_random.random(40, ["Cr", "Cu"], [3, 1])
