import os, sys, glob
import ROOT
from podio import root_io
import dd4hep as dd4hepModule
from ROOT import dd4hep

from math import sqrt

ROOT.gROOT.SetBatch(True)
#ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetTextSize(22)

bkg_file_folder = "CLD_SIM_nominal_configuration"
#bkg_file_folder = "CLD_SIM_low_energy_x_rays_EMZ_withAuger_deexcitationIgnoreCut_True"
plot_folder = bkg_file_folder + "_plots"
if not os.path.isdir(plot_folder):
    os.mkdir(plot_folder)

bkg_file_base_directory = "/eos/experiment/fcc/users/b/brfranco/background_files/guineaPig_andrea_June2024_v23_vtx000/"
bkg_file_pattern = os.path.join(bkg_file_base_directory, bkg_file_folder, "output_*.root")

file_list = glob.glob(bkg_file_pattern)
print(file_list)
print("Analyzing %d files from the following pattern %s"%(len(file_list), bkg_file_pattern))

def draw_histo(histo, drawOptions = "", path = plot_folder, logY = False, logX = False):
    name = histo.GetName()
    canvas = ROOT.TCanvas(name, name)
    if "colz" in drawOptions:
        histo.SetStats(0)
        canvas.SetRightMargin(0.15)
        if "map" in name:
            canvas.SetGrid()
        if "logz" in drawOptions:
            canvas.SetLogz()
            name += "_logz"
    histo.Draw(drawOptions)
    if logY:
        canvas.SetLogy()
        name += "_logy"
    if logX:
        canvas.SetLogx()
        name += "_logx"
    canvas.Print(os.path.join(path, name + ".png"))



#Declare TH1s
number_of_hit_vtx_barrel =  ROOT.TH1F(f"number_of_hit_vtx_barrel", f"number_of_hit_vtx_barrel", 50, 0, 150)
number_of_hit_vtx_barrel.SetTitle(f"Number of hit in the vertex barrel detector")
number_of_hit_vtx_barrel.GetXaxis().SetTitle("Number of hit in the vertex barrel detector")
number_of_hit_vtx_barrel.GetYaxis().SetTitle("Number of entries")

number_of_hit_above_1_kev_vtx_barrel =  ROOT.TH1F(f"number_of_hit_above_1_kev_vtx_barrel", f"number_of_hit_above_1_kev_vtx_barrel", 50, 0, 150)
number_of_hit_above_1_kev_vtx_barrel.SetTitle(f"Number of hit above one keV in the vertex barrel detector")
number_of_hit_above_1_kev_vtx_barrel.GetXaxis().SetTitle("Number of hit above 1 keV in the vertex barrel detector")
number_of_hit_above_1_kev_vtx_barrel.GetYaxis().SetTitle("Number of entries")

vtx_barrel_hit_energy_deposit =  ROOT.TH1F(f"vtx_barrel_hit_energy_deposit", f"vtx_barrel_hit_energy_deposit", 100, 0, 100)
vtx_barrel_hit_energy_deposit.SetTitle(f"Vertex barrel hit energy deposit")
vtx_barrel_hit_energy_deposit.GetXaxis().SetTitle("Vertex barrel hit energy deposit [keV]")
vtx_barrel_hit_energy_deposit.GetYaxis().SetTitle("Number of entries")

podio_reader = root_io.Reader(file_list)
for event in podio_reader.get("events"):
    vertex_barrel_hits = event.get("VertexBarrelCollection")
    number_of_hit_vtx_barrel.Fill(vertex_barrel_hits.size())
    number_of_hit_above_1_kev_vtx_barrel_counter = 0
    for vertex_barrel_hit in vertex_barrel_hits:
        vtx_barrel_hit_energy_deposit.Fill(vertex_barrel_hit.getEDep()*1000000)
        if vertex_barrel_hit.getEDep()*1000000 > 1:
            number_of_hit_above_1_kev_vtx_barrel_counter += 1
    number_of_hit_above_1_kev_vtx_barrel.Fill(number_of_hit_above_1_kev_vtx_barrel_counter)


draw_histo(number_of_hit_vtx_barrel)
draw_histo(number_of_hit_above_1_kev_vtx_barrel)
draw_histo(vtx_barrel_hit_energy_deposit)
