import os, sys
import ROOT
from podio import root_io
import dd4hep as dd4hepModule
from ROOT import dd4hep

from math import sqrt

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetTextSize(22)


#bkg_file_directory = "/eos/experiment/fcc/users/b/brfranco/background_files/guineaPig_andrea_Apr2024/sim_noboost_noBfield/"
#bkg_file_directory = "/eos/experiment/fcc/users/b/brfranco/background_files/guineaPig_andrea_Apr2024/sim_noboost_bfield/"
#bkg_file_directory = "/eos/experiment/fcc/users/b/brfranco/background_files/guineaPig_andrea_Apr2024/sim_boost_noBfield/"
#bkg_file_directory = "/eos/experiment/fcc/users/b/brfranco/background_files/guineaPig_andrea_Apr2024/sim_boost_bfield/"
#bkg_file_directory =  "/eos/experiment/fcc/users/b/brfranco/background_files/guineaPig_andrea_Apr2024/sim_boost_bfield_k4geo_version/"
bkg_file_directory =  "/eos/experiment/fcc/users/b/brfranco/background_files/guineaPig_andrea_Apr2024/sim_boost_bfield_k4geo_version_Geant4TrackerAction_edep0/"
file_name_template = "IDEA_01_v03_pairs_XX.root"
plot_folder = "Geant4TrackerAction_edep0_stepLength1m_800_ns"
number_of_iteration_on_bx_batches = 1
number_of_bx = 40
bunch_spacing = 20 # ns
if not os.path.isdir(plot_folder):
    os.mkdir(plot_folder)

def draw_histo(histo, drawOptions = "", path = plot_folder, logY = False):
    name = histo.GetName()
    canvas = ROOT.TCanvas(name, name)
    histo.Draw(drawOptions)
    if logY:
        canvas.SetLogy()
    canvas.Print(os.path.join(path, name + ".png"))

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
        total_number_of_layers += 1
        total_number_of_cells += n_cell_superlayer0 + sl * n_cell_increment
        n_cell_per_layer[str(n_layers_per_superlayer * sl + l)] = n_cell_superlayer0 + sl * n_cell_increment
print("total_number_of_cells: ", total_number_of_cells)
print("total_number_of_layers: ", total_number_of_layers)
print("n_cell_per_layer: ", n_cell_per_layer)

# Histogram definition
# what is the momentum spectrum of bkg particles (show tails, log scale)?
# on average, how many hit a bkg leave in the DC (provided that they lead to at least one hit): _DCHCollection_particle.index 
# th1 with x axis being the number of hit in the fired cells or a th1 with x axis being the layer indies, y axis being the number of hits received
# what is the energy of the particles actually hitting the DC
energy_of_particles_hitting_dch =  ROOT.TH1F(f"energy_of_particles_hitting_dch", f"energy_of_particles_hitting_dch", 50, 0, 3)
energy_of_particles_hitting_dch.SetTitle(f"Energy of bkg particles hitting the DCH, {number_of_iteration_on_bx_batches * number_of_bx} BXs")
energy_of_particles_hitting_dch.GetXaxis().SetTitle("Energy [GeV]")
energy_of_particles_hitting_dch.GetYaxis().SetTitle("Number of entries")
# what is the de/DX of the background particles/muons (from another rootifle)?
# what is the percentage of primaries/secondaries hitting the DC? Should be mostly secondaries
# what is the particle type hitting DC?
type_of_particles_hitting_dch =  ROOT.TH1F(f"type_of_particles_hitting_dch", f"type_of_particles_hitting_dch", 22, 0, 22)
type_of_particles_hitting_dch.SetTitle(f"Type of bkg particles hitting the DCH, {number_of_iteration_on_bx_batches * number_of_bx} BXs")
type_of_particles_hitting_dch.GetXaxis().SetTitle("PDG ID")
type_of_particles_hitting_dch.GetYaxis().SetTitle("Arbitrary units")

# how many particles hit the chamber per bx?

# where are the particles leading to hits in the DC coming from?
particle_hittingDC_origin_rz = ROOT.TH2F("particle_hittingDC_origin_rz", "particle_hittingDC_origin_rz", 100, 0, 2000, 100, 0, 200)
particle_hittingDC_origin_rz.SetTitle(f"Origin of bkg particles hitting the DCH ({number_of_iteration_on_bx_batches * number_of_bx} BXs)")
particle_hittingDC_origin_rz.GetXaxis().SetTitle("z [mm]")
particle_hittingDC_origin_rz.GetYaxis().SetTitle("r [mm]")


for bx_batch_index in range(0, number_of_iteration_on_bx_batches): # first loop to run X times on Y BX's
    print("Iteration number ", bx_batch_index)
    file_bx_index_list = [bx_batch_index * number_of_bx + i for i in range(1, number_of_bx + 1)]
    dict_cellID_nHits = {}
    total_number_of_hit_integrated = 0
    number_of_cell_with_multiple_hits = 0
    energy_deposit_R_z = ROOT.TH2F("energy_deposit_R_z", "energy_deposit_R_z", 50, 0, 2100, 50, 0, 2100)

    # Histogram declarations
    # what is the occupancy per radial layer?
    occupancy_per_layer =  ROOT.TH1F(f"occupancy_per_layer_bx_batch_index_{bx_batch_index}", f"occupancy_per_layer_bx_batch_index_{bx_batch_index}", total_number_of_layers, 0, total_number_of_layers)
    occupancy_per_layer.SetTitle(f"Occupancy per layer, {number_of_bx} BXs ({number_of_bx * bunch_spacing} ns)")
    occupancy_per_layer.GetXaxis().SetTitle("Radial layer index")
    occupancy_per_layer.GetYaxis().SetTitle("Channel occupancy [%]")

    # Where do the particles hit the DC?
    DC_simhit_position_rz = ROOT.TH2F(f"DC_simhit_position_rz_bx_batch_index_{bx_batch_index}", "DC_simhit_position_rz_bx_batch_index_{bx_batch_index}", 100, 0, 2000, 100, 0, 2000)
    DC_simhit_position_rz.SetTitle(f"DCH simHit position ({number_of_bx} BXs)")
    DC_simhit_position_rz.GetXaxis().SetTitle("z [mm]")
    DC_simhit_position_rz.GetYaxis().SetTitle("r [mm]")

    for file_bx_index in file_bx_index_list: # loop over the set of file for this batch of BX's
        input_file_path = os.path.join(bkg_file_directory, file_name_template.replace("XX", str(file_bx_index)))
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
                if layer >= n_layers_per_superlayer or superlayer >= n_superlayers:
                    print("Error: layer or super layer index out of range")
                    print(f"Layer: {layer} while max layer is {n_layers_per_superlayer - 1}. Superlayer: {superlayer} while max superlayer is {n_superlayers - 1}.")
                unique_layer_index = superlayer * n_layers_per_superlayer + layer
                nphi = decoder.get(cellID, "nphi")
                cellID_unique_identifier = "SL_" + str(superlayer)  + "_L_" + str(layer) + "_nphi_" + str(nphi) 
                if not cellID_unique_identifier in dict_cellID_nHits.keys(): # the cell was not fired yet
                    occupancy_per_layer.Fill(unique_layer_index)
                    dict_cellID_nHits[cellID_unique_identifier] = 1
                else:
                    if(dict_cellID_nHits[cellID_unique_identifier] == 1):
                        number_of_cell_with_multiple_hits += 1
                    dict_cellID_nHits[cellID_unique_identifier] += 1
                # where are the particles leading to hits in the DC coming from?
                particle = dc_hit.getParticle()
                particle_fourvector = ROOT.Math.LorentzVector('ROOT::Math::PxPyPzM4D<double>')(particle.getMomentum().x, particle.getMomentum().y, particle.getMomentum().z, particle.getMass())
                particle_origin = particle.getVertex()
                particle_hittingDC_origin_rz.Fill(abs(particle_origin.z), sqrt(particle_origin.x ** 2 + particle_origin.y ** 2))
                # Where do the particles hit the DC?
                DC_simhit_position_rz.Fill(abs(dc_hit.getPosition().z), sqrt(dc_hit.getPosition().x ** 2 + dc_hit.getPosition().y ** 2))
                # what is the particle type hitting DC?
                type_of_particles_hitting_dch.Fill(particle.getPDG())
                # what is the energy of the particles actually hitting the DC
                energy_of_particles_hitting_dch.Fill(particle_fourvector.E())

                # end of loop on sim hits

            # end of loop on events

        print("\t\tTotal number of bkg hit in the DC accumulating BXs: ", str(total_number_of_hit_integrated))
        print("\t\tTotal number of bkg hit in the DC from this BX: ", str(total_number_of_hit_thisbx))
        print("\t\tNumber of fired cells: ", str(len(dict_cellID_nHits.keys())))
        print("\t\tNumber of cells with more than one hit: ", str(number_of_cell_with_multiple_hits))
        print("\t\tPercentage of fired cells: ", str(len(dict_cellID_nHits.keys())/float(total_number_of_cells)))
    with open(os.path.join(plot_folder, "overall_occupancy.txt"), "w") as myfile:
        myfile.write(f"BX batch iteration {bx_batch_index}, integrating over {number_of_bx} BXs")
        myfile.write("\tTotal number of bkg hit in the DC accumulating BXs: " + str(total_number_of_hit_integrated))
        myfile.write("\tNumber of fired cells: " + str(len(dict_cellID_nHits.keys())))
        myfile.write("\tPercentage of fired cells: " + str(len(dict_cellID_nHits.keys())/float(total_number_of_cells)))
    print(os.path.join(plot_folder, "overall_occupancy.txt"), " written.")

    # Write the "per batch" histograms that have to be integrated over the BXs inside a single batch
    draw_histo(DC_simhit_position_rz, "colz")
    # Normalize the occupancy per layer th1
    for unique_layer_index in range(0, total_number_of_layers):
        raw_bin_content = occupancy_per_layer.GetBinContent(unique_layer_index + 1) # NB: we use the trick that the bin index here is the same as the layer index it corresponds to, just binIdx 0 is underflow
        occupancy_per_layer.SetBinContent(unique_layer_index + 1, raw_bin_content/float(n_cell_per_layer[str(unique_layer_index)])) # unique_layer_index and n_cell_per_layer key definitions coincide
    draw_histo(occupancy_per_layer)
    # end of loop on BX's for the given BX batch iteration

# end of loop on BX batches
type_of_particles_hitting_dch.Scale(1/type_of_particles_hitting_dch.Integral())
draw_histo(particle_hittingDC_origin_rz, "colz")
draw_histo(energy_of_particles_hitting_dch, logY = True)
draw_histo(type_of_particles_hitting_dch, logY = True)
