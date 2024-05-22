import os, sys
import ROOT
from podio import root_io
import dd4hep as dd4hepModule
from ROOT import dd4hep

from math import sqrt

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)

#bkg_file_directory = "/eos/experiment/fcc/users/b/brfranco/background_files/guineaPig_andrea_Apr2024/sim_noboost_noBfield/"
#bkg_file_directory = "/eos/experiment/fcc/users/b/brfranco/background_files/guineaPig_andrea_Apr2024/sim_noboost_bfield/"
#bkg_file_directory = "/eos/experiment/fcc/users/b/brfranco/background_files/guineaPig_andrea_Apr2024/sim_boost_noBfield/"
bkg_file_directory = "/eos/experiment/fcc/users/b/brfranco/background_files/guineaPig_andrea_Apr2024/sim_boost_bfield/"
file_name_template = "IDEA_01_v03_pairs_XX.root"
plot_folder = "baseline_400_ns"
if not os.path.isdir(plot_folder):
    os.mkdir(plot_folder)

# Drift chamber geometry parameters
n_layers_per_superlayer = 8
n_superlayers = 14
total_number_of_layers = 0
n_cell_superlayer0 = 192
n_cell_increment = 48
n_cell_per_layer = {}
total_number_of_cells = 0
for sl in range(0, n_superlayers):
    for l in range(0, n_layers_per_superlayer):
        total_number_of_layer += 1
        total_number_of_cells += n_cell_superlayer0 + sl * n_cell_increment
        n_cell_per_layer[str(n_layers_per_superlayer * sl + l)] = n_cell_superlayer0 + sl * n_cell_increment
print("total_number_of_cells: ", total_number_of_cells)
print("total_number_of_layers: ", total_number_of_layers)
print("n_cell_per_layer: ", n_cell_per_layer)

# Histogram definition
# where are the bkg particles coming from?
particle_origin_rz = ROOT.TH2F("particle_origin_rz", "particle_origin_rz", 500, 0, 2500, 500, 0, 2000)
primaries_origin_r = ROOT.TH1F("primaries_origin_r", "primaries_origin_r", 50, 0, 50)
# what is the momentum spectrum of bkg particles (show tails, log scale)?
# on average, how many hit a bkg leave in the DC (provided that they lead to at least one hit)
# th1 with x axis being the number of hit in the fired cells or a th1 with x axis being the layer indies, y axis being the number of hits received
# what is the momentum of the particles actually hitting the DC?
# what is the de/DX of the background particles/muons (from another rootifle)?
# what is the percentage of primaries/secondaries hitting the DC?
# what is the particle type hitting DC?

# where are the particles leading to hits in the DC coming from?
particle_hittingDC_origin_rz = ROOT.TH2F("particle_hittingDC_origin_rz", "particle_hittingDC_origin_rz", 100, 0, 2000, 100, 0, 200)
particle_hittingDC_origin_rz.GetXaxis().SetTitle("z [mm]")
particle_hittingDC_origin_rz.GetYaxis().SetTitle("r [mm]")

number_of_iteration = 1
number_of_bx = 40
for iteration in range(0, number_of_iteration): # first loop to run 5 times on 20 BX
    print("Iteration number ", iteration)
    index_list = [iteration * number_of_bx + i for i in range(1, number_of_bx + 1)]
    dict_cellID_nHits = {}
    total_number_of_hit_integrated = 0
    number_of_cell_with_multiple_hits = 0
    energy_deposit_R_z = ROOT.TH2F("energy_deposit_R_z", "energy_deposit_R_z", 50, 0, 2100, 50, 0, 2100)

    # Histogram declarations
    # what is the occupancy per radial layer?
    occupancy_per_layer =  ROOT.TH1I("occupancy_per_layer_iteration_" + str(iteration), "occupancy_per_layer_iteration_" + str(iteration), total_number_of_layers, 0, total_number_of_layers)

    for index in index_list: # loop over the set of file for this batch of BX's
        input_file_path = os.path.join(bkg_file_directory, file_name_template.replace("XX", str(index)))
        print("\tTreating: %s"%input_file_path)
        podio_reader = root_io.Reader(input_file_path)
        metadata = podio_reader.get("metadata")[0]
        cellid_encoding = metadata.get_parameter("DCHCollection__CellIDEncoding")
        decoder = dd4hep.BitFieldCoder(cellid_encoding)
        total_number_of_hit_thisbx = 0
        for event in podio_reader.get("events"):
            for dc_hit in event.get("DCHCollection"):
                total_number_of_hit_integrated += 1
                total_number_of_hit_thisbx += 1
                cellID = dc_hit.getCellID()
                layer = decoder.get(cellID, "layer")
                superlayer = decoder.get(cellID, "superlayer")
                unique_layer_index = superlayer * n_layers_per_superlayer + layer
                nphi = decoder.get(cellID, "nphi")
                cellID_unique_identifier = "SL_" + str(superlayer)  + "_L_" + str(layer) + "_nphi_" + str(nphi) 
                # check if the cell was fired already, if yes, just increment counter
                if not cellID_unique_identifier in dict_cellID_nHits.keys():
                    dict_cellID_nHits[cellID_unique_identifier] = 1
                else:
                    if(dict_cellID_nHits[cellID_unique_identifier] == 1):
                        number_of_cell_with_multiple_hits += 1
                    dict_cellID_nHits[cellID_unique_identifier] += 1
                # where are the particles leading to hits in the DC coming from?
                particle = dc_hit.getParticle()
                particle_origin = particle.getVertex()
                particle_hittingDC_origin_rz.Fill(abs(particle_origin.z), sqrt(particle_origin.x ** 2 + particle_origin.y ** 2))
        print("\t\tTotal number of bkg hit in the DC accumulating BXs: ", str(total_number_of_hit_integrated))
        print("\t\tTotal number of bkg hit in the DC from this BX: ", str(total_number_of_hit_thisbx))
        print("\t\tNumber of fired cells: ", str(len(dict_cellID_nHits.keys())))
        print("\t\tNumber of cells with more than one hit: ", str(number_of_cell_with_multiple_hits))
        print("\t\tPercentage of fired cells: ", str(len(dict_cellID_nHits.keys())/float(total_number_of_cells)))
    # end of loop on BX's

    # Fill the occupancy per layer (assuming a cell is fired even if it had more than one hit)
    for cellID in dict_cellID_nHits.keys():




def draw_histo(histo, path="./", drawOptions = ""):
    name = histo.GetName()
    canvas = ROOT.TCanvas(name, name)
    histo.Draw(drawOptions)
    canvas.Print(os.path.join(path, name + ".png"))

draw_histo(particle_hittingDC_origin_rz, plot_folder, "colz")
