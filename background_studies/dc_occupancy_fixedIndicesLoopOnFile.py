import os, sys
import ROOT
from podio import root_io
import dd4hep as dd4hepModule
from ROOT import dd4hep

from math import sqrt

ROOT.gROOT.SetBatch(True)
#ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetTextSize(22)


#bkg_file_directory = "/eos/experiment/fcc/users/b/brfranco/background_files/guineaPig_andrea_Apr2024/sim_noboost_noBfield/"
#bkg_file_directory = "/eos/experiment/fcc/users/b/brfranco/background_files/guineaPig_andrea_Apr2024/sim_noboost_bfield/"
#bkg_file_directory = "/eos/experiment/fcc/users/b/brfranco/background_files/guineaPig_andrea_Apr2024/sim_boost_noBfield/"
#bkg_file_directory = "/eos/experiment/fcc/users/b/brfranco/background_files/guineaPig_andrea_Apr2024/sim_boost_bfield/"
#bkg_file_directory =  "/eos/experiment/fcc/users/b/brfranco/background_files/guineaPig_andrea_Apr2024/sim_boost_bfield_k4geo_version/"
#bkg_file_directory =  "/eos/experiment/fcc/users/b/brfranco/background_files/guineaPig_andrea_Apr2024/sim_boost_bfield_k4geo_version_Geant4TrackerAction_edep0/"
# NEW beampie 
bkg_file_directory = "/eos/experiment/fcc/users/b/brfranco/background_files/guineaPig_andrea_June2024_v23/ddsimoutput/new_beampipe_default_ddsim_config/"
file_name_template = "IDEA_01_v03_pairs_XX.root"
plot_folder = "new_beampipe_default_ddsim_config_plots"
number_of_iteration_on_bx_batches = 2
number_of_bx = 40
bunch_spacing = 20 # ns
if not os.path.isdir(plot_folder):
    os.mkdir(plot_folder)

def draw_histo(histo, drawOptions = "", path = plot_folder, logY = False):
    name = histo.GetName()
    canvas = ROOT.TCanvas(name, name)
    if "colz" in drawOptions:
        histo.SetStats(0)
        canvas.SetRightMargin(0.15)
        if "map" in name:
            canvas.SetGrid()
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
max_n_cell_per_layer = n_cell_per_layer[str(total_number_of_layers - 1)]

# Histogram definition
# what is the energy of the particles actually hitting the DC
energy_of_particles_hitting_dch =  ROOT.TH1F(f"energy_of_particles_hitting_dch", f"energy_of_particles_hitting_dch", 50, 0, 0.1)
energy_of_particles_hitting_dch.SetTitle(f"Energy of bkg particles hitting the DCH ({number_of_iteration_on_bx_batches * number_of_bx} BXs)")
energy_of_particles_hitting_dch.GetXaxis().SetTitle("Energy [GeV]")
energy_of_particles_hitting_dch.GetYaxis().SetTitle("Number of entries")
# what is the energy of the particles actually hitting the DC when they come from photons
energy_of_particles_hitting_dch_with_photon_parent =  ROOT.TH1F(f"energy_of_particles_hitting_dch_with_photon_parent", f"energy_of_particles_hitting_dch_with_photon_parent", 50, 0, 0.1)
energy_of_particles_hitting_dch_with_photon_parent.SetTitle(f"Energy of photon originating particles hitting the DCH ({number_of_iteration_on_bx_batches * number_of_bx} BXs)")
energy_of_particles_hitting_dch_with_photon_parent.GetXaxis().SetTitle("Energy [GeV]")
energy_of_particles_hitting_dch_with_photon_parent.GetYaxis().SetTitle("Number of entries")
# what is the de/DX of the background particles/muons (from another rootfile)?
dedx_of_dch_hits =  ROOT.TH1F(f"dedx_of_dch_hits", f"dedx_of_dch_hits", 50, 0, 1)
dedx_of_dch_hits.SetTitle(f"dE/dx of DCH hits")
dedx_of_dch_hits.GetXaxis().SetTitle("dE/dx [keV/mm]")
dedx_of_dch_hits.GetYaxis().SetTitle("Number of entries")
# what is the particle pdgid hitting DC?
pdgid_of_particles_hitting_dch =  ROOT.TH1F(f"pdgid_of_particles_hitting_dch", f"pdgid_of_particles_hitting_dch", 46, -23, 23)
pdgid_of_particles_hitting_dch.SetTitle(f"PDG ID of bkg particles hitting the DCH")
pdgid_of_particles_hitting_dch.GetXaxis().SetTitle("PDG ID")
pdgid_of_particles_hitting_dch.GetYaxis().SetTitle("Arbitrary units")
# all particles leaving signals in DCH are electrons or positrons, is their parent a photon?
hasPhotonParent_of_particles_hitting_dch =  ROOT.TH1F(f"hasPhotonParent_of_particles_hitting_dch", f"hasPhotonParent_of_particles_hitting_dch", 2, 0, 2)
hasPhotonParent_of_particles_hitting_dch.SetTitle(f"Was the particle produced by a photon?")
hasPhotonParent_of_particles_hitting_dch.GetXaxis().SetTitle("Has photon parent?")
hasPhotonParent_of_particles_hitting_dch.GetYaxis().SetTitle("Arbitrary units")
# are the particles hitting dc primaries or secondaries?
isPrimary_of_particles_hitting_dch =  ROOT.TH1F(f"isPrimary_of_particles_hitting_dch", f"isPrimary_of_particles_hitting_dch", 2, 0, 2)
isPrimary_of_particles_hitting_dch.SetTitle(f"Is the particle a primary?")
isPrimary_of_particles_hitting_dch.GetXaxis().SetTitle("Is primary?")
isPrimary_of_particles_hitting_dch.GetYaxis().SetTitle("Arbitrary units")
# what is the vertex radius of the oldest parent (original primary particle)?
oldestParentVertexRadius_of_particles_hitting_dch =  ROOT.TH1F(f"oldestParentVertexRadius_of_particles_hitting_dch", f"oldestParentVertexRadius_of_particles_hitting_dch", 50, 0, 50)
oldestParentVertexRadius_of_particles_hitting_dch.SetTitle(f"Radius of oldest parent for particles hitting the DCH")
oldestParentVertexRadius_of_particles_hitting_dch.GetXaxis().SetTitle("Radius [mm]")
oldestParentVertexRadius_of_particles_hitting_dch.GetYaxis().SetTitle("Number of entries")
# where are the particles leading to hits in the DC coming from?
particle_hittingDC_origin_rz = ROOT.TH2F("particle_hittingDC_origin_rz", "particle_hittingDC_origin_rz", 100, 0, 2000, 100, 0, 200)
particle_hittingDC_origin_rz.SetTitle(f"Origin of bkg particles hitting the DCH ({number_of_iteration_on_bx_batches * number_of_bx} BXs)")
particle_hittingDC_origin_rz.GetXaxis().SetTitle("z [mm]")
particle_hittingDC_origin_rz.GetYaxis().SetTitle("r [mm]")
# how many cells are fired by a bkg particle?
n_cell_fired_of_particles_hitting_dch =  ROOT.TH1F(f"n_cell_fired_of_particles_hitting_dch", f"n_cell_fired_of_particles_hitting_dch", 50, 0, 50)
n_cell_fired_of_particles_hitting_dch.SetTitle(f"Number of cell fired by particles hitting the DCH")
n_cell_fired_of_particles_hitting_dch.GetXaxis().SetTitle("Number of cell fired by the particle")
n_cell_fired_of_particles_hitting_dch.GetYaxis().SetTitle("Number of entries")
# how many cells are fired by a bkg particle in log scale?
n_cell_fired_of_particles_hitting_dch_log =  ROOT.TH1F(f"n_cell_fired_of_particles_hitting_dch_log", f"n_cell_fired_of_particles_hitting_dch_log", 50, 0, 112)
n_cell_fired_of_particles_hitting_dch_log.SetTitle(f"Number of cell fired by particles hitting the DCH")
n_cell_fired_of_particles_hitting_dch_log.GetXaxis().SetTitle("Number of cell fired by the particle")
n_cell_fired_of_particles_hitting_dch_log.GetYaxis().SetTitle("Number of entries")

total_number_of_hit = 0
total_number_of_hit_comingFromPrimaryAfterBeamPipe = 0

string_for_overall_occupancy = ""

for bx_batch_index in range(0, number_of_iteration_on_bx_batches): # first loop to run X times on Y BX's
    print("Iteration number ", bx_batch_index)
    file_bx_index_list = [bx_batch_index * number_of_bx + i for i in range(1, number_of_bx + 1)]
    dict_cellID_nHits = {}
    total_number_of_hit_integrated_per_batch = 0
    number_of_cell_with_multiple_hits = 0

    # what fraction of particles do reach the drift chamber?
    total_number_of_particles_reaching_dch_per_bx_batch = 0
    total_number_of_particles_per_bx_batch = 0

    # Histogram declarations
    # what is the occupancy per radial layer?
    occupancy_per_layer =  ROOT.TH1F(f"occupancy_per_layer_bx_batch_index_{bx_batch_index}", f"occupancy_per_layer_bx_batch_index_{bx_batch_index}", total_number_of_layers, 0, total_number_of_layers)
    occupancy_per_layer.SetTitle(f"Occupancy per layer, {number_of_bx} BXs ({number_of_bx * bunch_spacing} ns)")
    occupancy_per_layer.GetXaxis().SetTitle("Radial layer index")
    occupancy_per_layer.GetYaxis().SetTitle("Channel occupancy [%]")

    # what is the energy deposited in cells, averaged over phi, per radial layer? #TODO
    cellEnergy_per_layer =  ROOT.TH1F(f"cellEnergy_per_layer_bx_batch_index_{bx_batch_index}", f"cellEnergy_per_layer_bx_batch_index_{bx_batch_index}", total_number_of_layers, 0, total_number_of_layers)
    cellEnergy_per_layer.SetTitle(f"Cell energy per layer, {number_of_bx} BXs ({number_of_bx * bunch_spacing} ns)")
    cellEnergy_per_layer.GetXaxis().SetTitle("Radial layer index")
    cellEnergy_per_layer.GetYaxis().SetTitle("Energy per cell [MeV]")

    # Where do the particles hit the DC?
    DC_simhit_position_rz = ROOT.TH2F(f"DC_simhit_position_rz_bx_batch_index_{bx_batch_index}", "DC_simhit_position_rz_bx_batch_index_{bx_batch_index}", 100, 0, 2000, 100, 0, 2000)
    DC_simhit_position_rz.SetTitle(f"DCH simHit position ({number_of_bx} BXs)")
    DC_simhit_position_rz.GetXaxis().SetTitle("z [mm]")
    DC_simhit_position_rz.GetYaxis().SetTitle("r [mm]")
    DC_simhit_position_rz.GetZaxis().SetTitle("Number of entries")

    # Map of the fired cells energies
    DC_fired_cell_map = ROOT.TH2F(f"DC_fired_cell_map_bx_batch_index_{bx_batch_index}", "DC_fired_cell_map_bx_batch_index_{bx_batch_index}", max_n_cell_per_layer, 0, max_n_cell_per_layer, total_number_of_layers, 0, total_number_of_layers) 
    DC_fired_cell_map.SetTitle(f"R-phi map of fired cells (energy in MeV on z axis) ({number_of_bx} BXs)")
    DC_fired_cell_map.GetXaxis().SetTitle("Cell phi index")
    DC_fired_cell_map.GetYaxis().SetTitle("Cell layer index")
    DC_fired_cell_map.GetZaxis().SetTitle("Energy [MeV]")

    for file_bx_index in file_bx_index_list: # loop over the set of file for this batch of BX's
        input_file_path = os.path.join(bkg_file_directory, file_name_template.replace("XX", str(file_bx_index)))
        print("\tTreating: %s"%input_file_path)
        podio_reader = root_io.Reader(input_file_path)
        metadata = podio_reader.get("metadata")[0]
        cellid_encoding = metadata.get_parameter("DCHCollection__CellIDEncoding")
        decoder = dd4hep.BitFieldCoder(cellid_encoding)
        for event in podio_reader.get("events"):
            total_number_of_hit_thisbx = 0
            seen_particle_ids = []
            dict_particle_n_fired_cell = {}
            dict_particle_fired_cell_id = {}
            total_number_of_particles_per_bx_batch += len(event.get("MCParticles"))
            for dc_hit in event.get("DCHCollection"):
                total_number_of_hit += 1
                total_number_of_hit_integrated_per_batch += 1
                total_number_of_hit_thisbx += 1
                cellID = dc_hit.getCellID()
                layer = decoder.get(cellID, "layer")
                superlayer = decoder.get(cellID, "superlayer")
                nphi = decoder.get(cellID, "nphi")
                particle = dc_hit.getParticle()
                if layer >= n_layers_per_superlayer or superlayer >= n_superlayers:
                    print("Error: layer or super layer index out of range")
                    print(f"Layer: {layer} while max layer is {n_layers_per_superlayer - 1}. Superlayer: {superlayer} while max superlayer is {n_superlayers - 1}.")
                unique_layer_index = superlayer * n_layers_per_superlayer + layer
                cellID_unique_identifier = "SL_" + str(superlayer)  + "_L_" + str(layer) + "_nphi_" + str(nphi) 
                # what is the occupancy?
                if not cellID_unique_identifier in dict_cellID_nHits.keys(): # the cell was not fired yet
                    occupancy_per_layer.Fill(unique_layer_index)
                    dict_cellID_nHits[cellID_unique_identifier] = 1
                else:
                    if(dict_cellID_nHits[cellID_unique_identifier] == 1):
                        number_of_cell_with_multiple_hits += 1
                    dict_cellID_nHits[cellID_unique_identifier] += 1
                # deal with the number of cell fired per particle #FIXME check that it does the proper thing
                particle_object_id = str(particle.getObjectID())
                if particle_object_id not in dict_particle_n_fired_cell.keys(): # the particle was not seen yet
                    dict_particle_n_fired_cell[particle_object_id] = 1
                    dict_particle_fired_cell_id[particle_object_id] = [cellID_unique_identifier]
                else: # the particle already fired cells
                    if not cellID_unique_identifier in dict_particle_fired_cell_id[particle_object_id]: # this cell was not yet fired by this particle
                        dict_particle_n_fired_cell[particle_object_id] += 1
                        dict_particle_fired_cell_id[particle_object_id].append(cellID_unique_identifier)
                # Where do the particles hit the DC?
                DC_simhit_position_rz.Fill(abs(dc_hit.getPosition().z), sqrt(dc_hit.getPosition().x ** 2 + dc_hit.getPosition().y ** 2))
                # Map of the fired cells energies
                DC_fired_cell_map.Fill(nphi, unique_layer_index, 1e+3*dc_hit.getEDep())
                # what is the de/DX of the background particles/muons (from another rootfile)?
                dedx_of_dch_hits.Fill(1e+6 * dc_hit.getEDep()/dc_hit.getPathLength())
                # fill the particle related TH1's (only once, a particle can lead to multiple DC hits)
                if particle.getObjectID() not in seen_particle_ids:
                    total_number_of_particles_reaching_dch_per_bx_batch += 1
                    # where are the particles leading to hits in the DC coming from?
                    particle_origin = particle.getVertex()
                    particle_hittingDC_origin_rz.Fill(abs(particle_origin.z), sqrt(particle_origin.x ** 2 + particle_origin.y ** 2))
                    # what is the energy of the particles actually hitting the DC
                    particle_fourvector = ROOT.Math.LorentzVector('ROOT::Math::PxPyPzM4D<double>')(particle.getMomentum().x, particle.getMomentum().y, particle.getMomentum().z, particle.getMass())
                    energy_of_particles_hitting_dch.Fill(particle_fourvector.E())
                    # what is the particle pdgid hitting DC?
                    pdgid_of_particles_hitting_dch.Fill(particle.getPDG())
                    # all particles leaving signals in DCH are electrons or positrons, is their parent a photon?
                    has_photon_parent = 0
                    for parent in particle.getParents():
                        if parent.getPDG() == 22:
                            has_photon_parent = 1
                            energy_of_particles_hitting_dch_with_photon_parent.Fill(particle_fourvector.E())
                    hasPhotonParent_of_particles_hitting_dch.Fill(has_photon_parent)

                    # are the particles hitting dc primaries or secondaries?
                    isPrimary_of_particles_hitting_dch.Fill(particle.getGeneratorStatus())
                    seen_particle_ids.append(particle.getObjectID()) # must be at the end

                # what is the vertex radius of the oldest parent (original primary particle)?
                # put this one outside of the seen_particle_ids condition to "weight" with the number of hit the particle lead to
                is_oldest_parent = False
                current_parent = particle
                while not is_oldest_parent:
                    if current_parent.parents_size() != 0:
                        current_parent = current_parent.getParents(0)
                    else:
                        is_oldest_parent = True
                # what is the vertex radius of the oldest parent (original primary particle)?
                parent_origin = current_parent.getVertex()
                parent_vertex_radius = sqrt(parent_origin.x ** 2 + parent_origin.y ** 2)
                oldestParentVertexRadius_of_particles_hitting_dch.Fill(parent_vertex_radius)
                if parent_vertex_radius > 10:
                    total_number_of_hit_comingFromPrimaryAfterBeamPipe += 1

                # end of loop on sim hits
            for particleKey in dict_particle_n_fired_cell.keys():
                n_cell_fired_of_particles_hitting_dch.Fill(dict_particle_n_fired_cell[particleKey])
                n_cell_fired_of_particles_hitting_dch_log.Fill(dict_particle_n_fired_cell[particleKey])

            print("\t\t\tTotal number of bkg hit in the DC from this BX: ", str(total_number_of_hit_thisbx))
            # end of loop on events

        print("\t\tTotal number of bkg hit in the DC accumulating BXs: ", str(total_number_of_hit_integrated_per_batch))
        print("\t\tNumber of fired cells: ", str(len(dict_cellID_nHits.keys())))
        print("\t\tNumber of cells with more than one hit: ", str(number_of_cell_with_multiple_hits))
        print("\t\tPercentage of fired cells: ", str(100 * len(dict_cellID_nHits.keys())/float(total_number_of_cells)), " %")
        print("\t\tTotal number of particles so far: ", total_number_of_particles_per_bx_batch)
        print("\t\tTotal number of particles that have hit the DC so far: ", total_number_of_particles_reaching_dch_per_bx_batch)
        print("\t\tPercentage of particles reaching the DC so far: ", 100 * total_number_of_particles_reaching_dch_per_bx_batch / float(total_number_of_particles_per_bx_batch), " %")

    string_for_overall_occupancy += f"BX batch iteration {bx_batch_index}, integrating over {number_of_bx} BXs\n"
    string_for_overall_occupancy += f"\tTotal number of bkg hit in the DC accumulating BXs: {total_number_of_hit_integrated_per_batch}\n"
    string_for_overall_occupancy += f"\tNumber of fired cells: {len(dict_cellID_nHits.keys())}\n"
    string_for_overall_occupancy += f"\tPercentage of fired cells: {len(dict_cellID_nHits.keys())/float(total_number_of_cells)}\n"

    # Write the "per batch" histograms that have to be integrated over the BXs inside a single batch
    draw_histo(DC_simhit_position_rz, "colz")
    draw_histo(DC_fired_cell_map, "colz")
    # Normalize the occupancy per layer th1
    for unique_layer_index in range(0, total_number_of_layers):
        raw_bin_content = occupancy_per_layer.GetBinContent(unique_layer_index + 1) # NB: we use the trick that the bin index here is the same as the layer index it corresponds to, just binIdx 0 is underflow
        occupancy_per_layer.SetBinContent(unique_layer_index + 1, 100 * raw_bin_content/float(n_cell_per_layer[str(unique_layer_index)])) # unique_layer_index and n_cell_per_layer key definitions coincide
    draw_histo(occupancy_per_layer)
    # end of loop on BX's for the given BX batch iteration

print(f"Percentage of hit coming from primary particles with vertex radius after beampipe: {total_number_of_hit_comingFromPrimaryAfterBeamPipe / float(total_number_of_hit)}")
# end of loop on BX batches
draw_histo(particle_hittingDC_origin_rz, "colz")
draw_histo(dedx_of_dch_hits)
draw_histo(energy_of_particles_hitting_dch, logY = True)
draw_histo(energy_of_particles_hitting_dch_with_photon_parent, logY = True)
pdgid_of_particles_hitting_dch.Scale(1/pdgid_of_particles_hitting_dch.Integral())
draw_histo(pdgid_of_particles_hitting_dch)
draw_histo(hasPhotonParent_of_particles_hitting_dch)
isPrimary_of_particles_hitting_dch.Scale(1/isPrimary_of_particles_hitting_dch.Integral())
draw_histo(isPrimary_of_particles_hitting_dch)
draw_histo(oldestParentVertexRadius_of_particles_hitting_dch)
draw_histo(n_cell_fired_of_particles_hitting_dch)
draw_histo(n_cell_fired_of_particles_hitting_dch_log, logY = True)

with open(os.path.join(plot_folder, "overall_occupancy.txt"), "w") as myfile:
    myfile.write(string_for_overall_occupancy)
print(os.path.join(plot_folder, "overall_occupancy.txt"), " written.")
