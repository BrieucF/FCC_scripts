import os, sys, glob
import ROOT
from podio import root_io
import dd4hep as dd4hepModule
from ROOT import dd4hep

from math import sqrt

ROOT.gROOT.SetBatch(True)
#ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetTextSize(22)

nominal_bkg_file_folder = "CLD_SIM_nominal_configuration"
#detailed_bkg_file_folder = "CLD_SIM_low_energy_x_rays_EMZ_withAuger_deexcitationIgnoreCut_True"
detailed_bkg_file_folder = "CLD_SIM_low_energy_x_rays_EMZ_witouthAuger"
plot_folder = detailed_bkg_file_folder + "_ratio_plots"
if not os.path.isdir(plot_folder):
    os.mkdir(plot_folder)


bkg_file_base_directory = "/eos/experiment/fcc/users/b/brfranco/background_files/guineaPig_andrea_June2024_v23_vtx000/"

file_name_list = os.listdir(os.path.join(bkg_file_base_directory, detailed_bkg_file_folder))

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
number_of_hit_vtx_barrel_ratio =  ROOT.TH1F(f"number_of_hit_vtx_barrel_ratio", f"number_of_hit_vtx_barrel_ratio", 50, 0, 3)
number_of_hit_vtx_barrel_ratio.SetTitle(f"Ratio of the number of hit in the vertex barrel detector: detailed simulation / nominal simulation")
number_of_hit_vtx_barrel_ratio.GetXaxis().SetTitle("NhitDetailed / NhitNominal")
number_of_hit_vtx_barrel_ratio.GetYaxis().SetTitle("Number of entries")

number_of_hit_above_1_kev_vtx_barrel_ratio =  ROOT.TH1F(f"number_of_hit_above_1_kev_vtx_barrel_ratio", f"number_of_hit_above_1_kev_vtx_barrel_ratio", 50, 0, 3)
number_of_hit_above_1_kev_vtx_barrel_ratio.SetTitle(f"Ratio of the number of hit > 1 keV in the vertex barrel detector: detailed simulation / nominal simulation")
number_of_hit_above_1_kev_vtx_barrel_ratio.GetXaxis().SetTitle("NhitAbove1keVDetailed / NhitAbove1keVNominal")
number_of_hit_above_1_kev_vtx_barrel_ratio.GetYaxis().SetTitle("Number of entries")

for file_name in file_name_list:
    nominal_file = os.path.join(bkg_file_base_directory, nominal_bkg_file_folder, file_name)
    nominal_podio_reader = root_io.Reader(nominal_file)
    nominal_event = nominal_podio_reader.get("events")[0]
    nominal_vertex_barrel_hits = nominal_event.get("VertexBarrelCollection")
    nominal_vertex_barrel_n_hits = nominal_vertex_barrel_hits.size()
    nominal_number_of_hit_above_1_kev_vtx_barrel_counter = 0
    for vertex_barrel_hit in nominal_vertex_barrel_hits:
        if vertex_barrel_hit.getEDep()*1000000 > 1:
            nominal_number_of_hit_above_1_kev_vtx_barrel_counter += 1
    del nominal_podio_reader

    detailed_file = os.path.join(bkg_file_base_directory, detailed_bkg_file_folder, file_name)
    detailed_podio_reader = root_io.Reader(detailed_file)
    detailed_event = detailed_podio_reader.get("events")[0]
    detailed_vertex_barrel_hits = detailed_event.get("VertexBarrelCollection")
    detailed_vertex_barrel_n_hits = detailed_vertex_barrel_hits.size()
    detailed_number_of_hit_above_1_kev_vtx_barrel_counter = 0
    for vertex_barrel_hit in detailed_vertex_barrel_hits:
        if vertex_barrel_hit.getEDep()*1000000 > 1:
            detailed_number_of_hit_above_1_kev_vtx_barrel_counter += 1
    del detailed_podio_reader

    number_of_hit_vtx_barrel_ratio.Fill(detailed_vertex_barrel_n_hits / nominal_vertex_barrel_n_hits)
    number_of_hit_above_1_kev_vtx_barrel_ratio.Fill(detailed_number_of_hit_above_1_kev_vtx_barrel_counter / nominal_number_of_hit_above_1_kev_vtx_barrel_counter)

draw_histo(number_of_hit_vtx_barrel_ratio)
draw_histo(number_of_hit_above_1_kev_vtx_barrel_ratio)

