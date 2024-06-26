import os
import glob

input_file_folder_template = "/eos/experiment/fcc/users/b/brfranco/background_files/guineaPig_andrea_June2024_v23/data*" #data2860/pairs.pairs"
output_file_doler = "/eos/experiment/fcc/users/b/brfranco/background_files/guineaPig_andrea_June2024_v23_vtx000/"

for data_folder_path in glob.glob(input_file_folder_template):
    bx_id = data_folder_path.split("/")[-1].replace("data", "")
    print(bx_id)
    input_file = os.path.join(data_folder_path, "pairs.pairs")
    if not os.path.exists(input_file):
        print(f"Skipping {bx_id}, file not produced previously")
        continue
    new_file_content = ""
    with open(input_file) as myfile:
        for line in myfile.readlines():
            entries = line.split(" ")
            # set the vtx to 000
            entries[4] = 0
            entries[5] = 0
            entries[6] = 0
            for entry in entries:
                new_file_content += str(entry) + " "
        new_file_content += "\n"

    output_file_folder = os.path.join(output_file_doler, f"data{bx_id}")
    if not os.path.isdir(output_file_folder):
        os.mkdir(output_file_folder)
    with open(os.path.join(output_file_folder, "pairs.pairs"), "w") as myoutputfile:
        myoutputfile.write(new_file_content)

    
