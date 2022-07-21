# -*- coding: utf-8 -*-

#  This script was made by CCD
#  AIRSS_v3.0_20220525
#  Condensed Matter Physics Group

#載入函數
from ase import Atom,Atoms
from ase.io.vasp import write_vasp,read_vasp
from ase.build import sort
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
		vec = positions1 - positions2						#兩位置向量差
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
			rs = Atoms(cell = ini_cell)

			#加入原子與位置
			for elements in append_list:
				condition = 0						#重置while True
				times = 0

				#如果rs是空的，在原點添加元素
				if len(rs.get_scaled_positions()) == 0:
					rs.append(elements)
					#print("first, ok")
					#print("now rs", rs, rs.get_scaled_positions())
					continue					#跳過後續進入下一個for迴圈

				#生成亂數並判斷距離
				while condition == 0:					#無限迴圈
					random = np.random.rand(1, 3)			#生成[1*3]亂數(位置)
					max_dis_times = 0				#小於最大距離次數
					times +=1
					#print(elements, times, "times random")
							
					#判斷距離
					for r_pos in range(len(rs.get_scaled_positions())):
						distance = Atoms_random.get_distance(random, rs.get_scaled_positions()[r_pos], rs.get_cell())
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
								atom2 = Atoms(elements, scaled_positions = random, cell = rs.get_cell())	#添加元素Atoms
								rs += atom2									#加到rs(原來的Atoms)
								rs.center(vacuum = add_vac)							#將cluster移到中心並加上真空層
								rs = sort(rs)									#sort symbol
								condition += 1									#使while False，打斷迴圈
								#print("ok, < max", max_dis_times, "times")

							else:				#最近原子不在最大距離內
								#print("retry, < max", max_dis_times, "times")
								break			#打斷for迴圈，回到while迴圈重新生成亂數
				
			#確認距離
			#Atoms_random.check_distance(rs.get_scaled_positions(), rs.get_cell())

			#調整真空層
			rs.center(vacuum = final_vac)

			#輸出vasp檔
			subprocess.run("mkdir " + str(rs.symbols) + "_" + "{:02d}".format(number + 1), shell = True)	#創建資料夾
			write_vasp(str(rs.symbols) + "_" + "{:02d}".format(number + 1) + "/POSCAR", rs, direct=True, sort=True, vasp5=True)
			print(number + 1, "complete")

<<<<<<< HEAD
Atoms_random.random(20, ["Cu", "Ni", "Si", "Cr"], [10, 10, 10, 10])
=======
Atoms_random.random(20, ["Cu", "Si"], [20, 20])
>>>>>>> AIRSS_V1.0
