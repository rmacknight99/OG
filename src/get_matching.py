import sys
import os



ID_to_match = sys.argv[1]

# Ignore for now
if '5x' in ID_to_match:
	sys.exit(0)

ID_to_match = ID_to_match.split('_')


is_list = False 
if len(ID_to_match) > 1:
	ID_to_match_1, ID_to_match_2 = int(ID_to_match[0]), int(ID_to_match[1])
	is_list = True 
else:
	ID_to_match = int(ID_to_match[0])


directory = os.fsencode('data/excited_state_monomer_spectra_csvs')
    
IDs = [file[:-13].decode("utf-8") for file in os.listdir(directory)]


count = 0 
matching = False
for ID in IDs:
	nums = ID.split('_')
	if is_list and len(nums) > 1:
		ID_1, ID_2 = int(nums[0]), int(nums[1])
		if (ID_1, ID_2) == (ID_to_match_1, ID_to_match_2) or (ID_2, ID_1) == (ID_to_match_1, ID_to_match_2):
			ID_to_write = ID
			matching = True 
			break
	elif not is_list and len(nums) == 1:
		ID_1 = nums[0]
		if int(ID_1) == int(ID_to_match):
			ID_to_write = str(ID_1)
			matching = True 
			break

if matching:
	with open('monomer_matches.csv', 'a') as f:
		f.write(sys.argv[1] + ',' + ID_to_write + '\n')




	
