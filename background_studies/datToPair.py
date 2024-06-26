import os

input_file_path = "/eos/experiment/fcc/users/b/brfranco/background_files/guineaPig_andrea_June2024_v23/"
for folder in os.listdir(input_file_path):
    input_filename = os.path.join(input_file_path, folder, "pairs.dat")
    if os.path.exists(input_filename):
        os.rename(input_filename, input_filename.replace(".dat", ".pairs"))

