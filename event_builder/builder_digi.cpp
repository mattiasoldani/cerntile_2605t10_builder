// DIGI-only event builder for runs without a FERS counterpart
// by M. Soldani, 2026 - developed with OpenAI Codex (GPT-5.5)

/* compilation:
g++ -std=c++17 builder_digi.cpp $(root-config --cflags --libs) -o builder_digi
*/

// usage: ./builder_digi <digi_root_file_or_dir> <digi_run_id> <output_dir> [recreate_outputs=1]
/* --> examples:
./builder_digi  /home/msoldani/25-27_cern/26_05_bt_t10/data_prep/data/TEST_root/digi_root/splitted  1778243915  ../data/TEST_root/global_root/splitted  1
./builder_digi  /eos/experiment/newtile/beamtests/26_05_t10/digi_root/splitted  1778243915  /eos/experiment/newtile/beamtests/26_05_t10/global_root/splitted  1
*/

#include <TDirectory.h>
#include <TFile.h>
#include <TKey.h>
#include <TLeaf.h>
#include <TObject.h>
#include <TString.h>
#include <TTimeStamp.h>
#include <TTree.h>

#include <algorithm>
#include <filesystem>
#include <iostream>
#include <stdexcept>
#include <string>
#include <system_error>
#include <vector>

#define NWFSAMPLES 1030
#define NFERSBDS 2
#define NFERSCHS1BD 64

namespace fs = std::filesystem;

namespace {

struct FersBuffers {
    Double_t TStamp_us[NFERSBDS] = {};
    Int_t Num_Hits[NFERSBDS] = {};
    Double_t dTRef[NFERSBDS] = {};
    ULong64_t Trg_Id[NFERSBDS] = {};
    ULong64_t ch_mask[NFERSBDS] = {};
    Short_t data_type[NFERSBDS][NFERSCHS1BD] = {};
    Int_t PHA_LG[NFERSBDS][NFERSCHS1BD] = {};
    Int_t PHA_HG[NFERSBDS][NFERSCHS1BD] = {};
    Float_t ToA[NFERSBDS][NFERSCHS1BD] = {};
    Float_t ToT[NFERSBDS][NFERSCHS1BD] = {};
    Int_t is_digi = 0;
    Int_t is_fers[NFERSBDS] = {};
};

struct DigiBuffers {
    Double_t run = 0;
    Double_t event0 = -1;
    Double_t trigger_ts = 0;
    Double_t event_id = -1;
    Long64_t trigger_ts_long = 0;
    Long64_t event_id_long = -1;
    bool trigger_ts_is_long = false;
    bool event_id_is_long = false;
    Double_t ro_time = 0;
    Double_t wave0[NWFSAMPLES] = {};
    Double_t wave1[NWFSAMPLES] = {};
    Double_t wave2[NWFSAMPLES] = {};
    Double_t wave3[NWFSAMPLES] = {};
    Double_t wave4[NWFSAMPLES] = {};
    Double_t wave5[NWFSAMPLES] = {};
    Double_t wave6[NWFSAMPLES] = {};
    Double_t wave7[NWFSAMPLES] = {};
    Int_t is_digi = 0;
    Int_t is_fers[NFERSBDS] = {};
};

struct InfoBuffers {
    UInt_t n_boards = 0;
    UInt_t n_channels = 0;
    UInt_t max_hits = 0;
    UShort_t board_mod = 0;
    TString file_format;
    TString janus_rel;
    TString acq_mode;
    UShort_t run = 0;
    ULong64_t digi_run = 0;
    ULong64_t time_epoch = 0;
    TTimeStamp time_UTC;
    UShort_t e_Nbins = 0;
    Double_t time_LSB_ns = 0;
    TString time_unit;
};

void collectTrees(TDirectory* dir, std::vector<TTree*>& trees)
{
    TIter next(dir->GetListOfKeys());
    while (auto* key = dynamic_cast<TKey*>(next())) {
        TObject* obj = key->ReadObj();
        if (!obj) continue;

        if (obj->InheritsFrom(TTree::Class())) {
            trees.push_back(dynamic_cast<TTree*>(obj));
            continue;
        }

        if (obj->InheritsFrom(TDirectory::Class())) {
            collectTrees(dynamic_cast<TDirectory*>(obj), trees);
        }
        delete obj;
    }
}

TTree* firstTree(TFile& file)
{
    std::vector<TTree*> trees;
    collectTrees(&file, trees);
    if (trees.empty()) {
        throw std::runtime_error(std::string("No TTree found in ") + file.GetName());
    }
    if (trees.size() > 1) {
        std::cerr << "Warning: found " << trees.size() << " TTrees in " << file.GetName()
                  << "; using '" << trees.front()->GetName() << "'.\n";
    }
    return trees.front();
}

void setDigiBranchAddresses(TTree& digi_tree, DigiBuffers& buffers)
{
    const std::string trigger_type = digi_tree.GetLeaf("trigger_ts")->GetTypeName();
    const std::string event_id_type = digi_tree.GetLeaf("event_id")->GetTypeName();
    buffers.trigger_ts_is_long = (trigger_type == "Long64_t");
    buffers.event_id_is_long = (event_id_type == "Long64_t");

    digi_tree.SetBranchAddress("run", &buffers.run);
    digi_tree.SetBranchAddress("event0", &buffers.event0);
    if (buffers.trigger_ts_is_long) {
        digi_tree.SetBranchAddress("trigger_ts", &buffers.trigger_ts_long);
    } else {
        digi_tree.SetBranchAddress("trigger_ts", &buffers.trigger_ts);
    }
    if (buffers.event_id_is_long) {
        digi_tree.SetBranchAddress("event_id", &buffers.event_id_long);
    } else {
        digi_tree.SetBranchAddress("event_id", &buffers.event_id);
    }
    digi_tree.SetBranchAddress("ro_time", &buffers.ro_time);
    digi_tree.SetBranchAddress("wave0", buffers.wave0);
    digi_tree.SetBranchAddress("wave1", buffers.wave1);
    digi_tree.SetBranchAddress("wave2", buffers.wave2);
    digi_tree.SetBranchAddress("wave3", buffers.wave3);
    digi_tree.SetBranchAddress("wave4", buffers.wave4);
    digi_tree.SetBranchAddress("wave5", buffers.wave5);
    digi_tree.SetBranchAddress("wave6", buffers.wave6);
    digi_tree.SetBranchAddress("wave7", buffers.wave7);
}

TTree* createDigiOutputTree(DigiBuffers& buffers)
{
    auto* tree = new TTree("digi", "All digi events");
    const std::string board_dim = "[" + std::to_string(NFERSBDS) + "]";
    const std::string wave_dim = "[" + std::to_string(NWFSAMPLES) + "]";
    tree->Branch("is_digi", &buffers.is_digi, "is_digi/I");
    tree->Branch("is_fers", buffers.is_fers, ("is_fers" + board_dim + "/I").c_str());
    tree->Branch("event0", &buffers.event0, "event0/D");
    tree->Branch("trigger_ts", &buffers.trigger_ts, "trigger_ts/D");
    tree->Branch("event_id", &buffers.event_id, "event_id/D");
    tree->Branch("ro_time", &buffers.ro_time, "ro_time/D");
    tree->Branch("wave0", buffers.wave0, ("wave0" + wave_dim + "/D").c_str());
    tree->Branch("wave1", buffers.wave1, ("wave1" + wave_dim + "/D").c_str());
    tree->Branch("wave2", buffers.wave2, ("wave2" + wave_dim + "/D").c_str());
    tree->Branch("wave3", buffers.wave3, ("wave3" + wave_dim + "/D").c_str());
    tree->Branch("wave4", buffers.wave4, ("wave4" + wave_dim + "/D").c_str());
    tree->Branch("wave5", buffers.wave5, ("wave5" + wave_dim + "/D").c_str());
    tree->Branch("wave6", buffers.wave6, ("wave6" + wave_dim + "/D").c_str());
    tree->Branch("wave7", buffers.wave7, ("wave7" + wave_dim + "/D").c_str());
    return tree;
}

TTree* createFersOutputTree(FersBuffers& buffers)
{
    auto* tree = new TTree("fers", "All merged FERS events");
    const std::string board_dim = "[" + std::to_string(NFERSBDS) + "]";
    const std::string board_channel_dims =
        "[" + std::to_string(NFERSBDS) + "][" + std::to_string(NFERSCHS1BD) + "]";
    tree->Branch("is_digi", &buffers.is_digi, "is_digi/I");
    tree->Branch("is_fers", buffers.is_fers, ("is_fers" + board_dim + "/I").c_str());
    tree->Branch("TStamp_us", buffers.TStamp_us, ("TStamp_us" + board_dim + "/D").c_str());
    tree->Branch("Num_Hits", buffers.Num_Hits, ("Num_Hits" + board_dim + "/I").c_str());
    tree->Branch("dTRef", buffers.dTRef, ("dTRef" + board_dim + "/D").c_str());
    tree->Branch("Trg_Id", buffers.Trg_Id, ("Trg_Id" + board_dim + "/l").c_str());
    tree->Branch("ch_mask", buffers.ch_mask, ("ch_mask" + board_dim + "/l").c_str());
    tree->Branch("data_type", buffers.data_type, ("data_type" + board_channel_dims + "/S").c_str());
    tree->Branch("PHA_LG", buffers.PHA_LG, ("PHA_LG" + board_channel_dims + "/I").c_str());
    tree->Branch("PHA_HG", buffers.PHA_HG, ("PHA_HG" + board_channel_dims + "/I").c_str());
    tree->Branch("ToA", buffers.ToA, ("ToA" + board_channel_dims + "/F").c_str());
    tree->Branch("ToT", buffers.ToT, ("ToT" + board_channel_dims + "/F").c_str());
    return tree;
}

void writeInfoTree(ULong64_t digi_run)
{
    InfoBuffers buffers;
    buffers.digi_run = digi_run;

    TTree info("info", "Event-independent information");
    info.Branch("DIGI_run", &buffers.digi_run, "DIGI_run/l");
    info.Branch("FERS_n_boards", &buffers.n_boards, "FERS_n_boards/i");
    info.Branch("FERS_n_channels", &buffers.n_channels, "FERS_n_channels/i");
    info.Branch("FERS_max_hits", &buffers.max_hits, "FERS_max_hits/i");
    info.Branch("FERS_board_mod", &buffers.board_mod, "FERS_board_mod/s");
    info.Branch("FERS_file_format", &buffers.file_format);
    info.Branch("FERS_janus_rel", &buffers.janus_rel);
    info.Branch("FERS_acq_mode", &buffers.acq_mode);
    info.Branch("FERS_run", &buffers.run, "FERS_run/s");
    info.Branch("FERS_time_epoch", &buffers.time_epoch, "FERS_time_epoch/l");
    info.Branch("FERS_time_UTC", &buffers.time_UTC, 32000, 0);
    info.Branch("FERS_e_Nbins", &buffers.e_Nbins, "FERS_e_Nbins/s");
    info.Branch("FERS_time_LSB_ns", &buffers.time_LSB_ns, "FERS_time_LSB_ns/D");
    info.Branch("FERS_time_unit", &buffers.time_unit);
    info.Fill();
    info.Write();
}

void copyDigi(const DigiBuffers& input, DigiBuffers& output)
{
    output = input;
    if (input.trigger_ts_is_long) output.trigger_ts = input.trigger_ts_long;
    if (input.event_id_is_long) output.event_id = input.event_id_long;
    output.is_digi = 1;
    for (int board = 0; board < NFERSBDS; ++board) {
        output.is_fers[board] = 0;
    }
}

void resetFers(FersBuffers& buffers)
{
    buffers.is_digi = 0;
    for (int board = 0; board < NFERSBDS; ++board) {
        buffers.TStamp_us[board] = 0;
        buffers.Num_Hits[board] = 0;
        buffers.dTRef[board] = 0;
        buffers.Trg_Id[board] = 0;
        buffers.ch_mask[board] = 0;
        buffers.is_fers[board] = 0;
        for (int ch = 0; ch < NFERSCHS1BD; ++ch) {
            buffers.data_type[board][ch] = 0;
            buffers.PHA_LG[board][ch] = -1;
            buffers.PHA_HG[board][ch] = -1;
            buffers.ToA[board][ch] = -1;
            buffers.ToT[board][ch] = -1;
        }
    }
}

bool belongsToRun(const fs::path& path, const std::string& digi_run_id)
{
    return path.filename().string().rfind(digi_run_id + "_", 0) == 0;
}

std::vector<std::string> expandDigiInputs(const std::string& digi_path, const std::string& digi_run_id)
{
    std::vector<std::string> files;
    if (fs::is_directory(digi_path)) {
        for (const auto& entry : fs::directory_iterator(digi_path)) {
            if (!entry.is_regular_file()) continue;
            if (entry.path().extension() != ".root") continue;
            if (!belongsToRun(entry.path(), digi_run_id)) continue;
            files.push_back(entry.path().string());
        }
    } else {
        files.push_back(digi_path);
    }

    std::sort(files.begin(), files.end());
    return files;
}

std::string digiSplitId(const std::string& digi_file_path, const std::string& digi_run_id)
{
    const std::string stem = fs::path(digi_file_path).stem().string();
    const std::string prefix = digi_run_id + "_";
    if (stem.rfind(prefix, 0) == 0) return stem.substr(prefix.size());
    return stem;
}

std::string outputFilePath(
    const std::string& output_path,
    const std::string& digi_run_id,
    const std::string& digi_file_path)
{
    const std::string name = "XXXX_" + digi_run_id + "_"
        + digiSplitId(digi_file_path, digi_run_id) + ".root";
    return (fs::path(output_path) / name).string();
}

bool parseRecreateOutputs(int argc, char** argv)
{
    if (argc == 4) return true;

    const std::string value = argv[4];
    if (value == "0") return false;
    if (value == "1") return true;
    throw std::runtime_error("recreate_outputs must be 0 or 1");
}

Long64_t buildOneFile(
    const std::string& digi_path,
    const std::string& digi_run_id,
    const std::string& output_path,
    bool recreate_outputs)
{
    const std::string out_path = outputFilePath(output_path, digi_run_id, digi_path);
    if (fs::exists(out_path)) {
        if (!recreate_outputs) {
            std::cout << "Keeping existing " << out_path << "; skipping.\n";
            return 0;
        } else {
            std::cout << "Recreating " << out_path << "...\n";
        }

        std::error_code ec;
        if (!fs::remove(out_path, ec) || ec) {
            throw std::runtime_error("Cannot remove existing output file: " + out_path);
        }
    } else {
        std::cout << "Creating (1st time) " << out_path << "...\n";
    }

    TFile digi_file(digi_path.c_str(), "READ");
    if (digi_file.IsZombie()) throw std::runtime_error("Cannot open digi ROOT file: " + digi_path);
    TTree* digi_tree = firstTree(digi_file);

    TFile out_file(out_path.c_str(), "CREATE");
    if (out_file.IsZombie()) throw std::runtime_error("Cannot create output ROOT file: " + out_path);

    DigiBuffers digi_input;
    DigiBuffers digi_output;
    FersBuffers fers_output;
    setDigiBranchAddresses(*digi_tree, digi_input);
    auto* out_digi = createDigiOutputTree(digi_output);
    auto* out_fers = createFersOutputTree(fers_output);

    const Long64_t entries = digi_tree->GetEntries();
    Long64_t written = 0;
    for (Long64_t entry = 0; entry < entries; ++entry) {
        digi_tree->GetEntry(entry);
        copyDigi(digi_input, digi_output);
        resetFers(fers_output);
        fers_output.is_digi = digi_output.is_digi;

        out_digi->Fill();
        out_fers->Fill();
        ++written;

        const Long64_t processed = entry + 1;
        if (processed % 10000 == 0 || processed == entries) {
            std::cout << "Processed " << processed << "/" << entries
                      << " input DIGI entries; in total: wrote " << written << ".\n";
        }
    }

    out_file.cd();
    out_digi->Write();
    out_fers->Write();
    writeInfoTree(std::stoull(digi_run_id));
    out_file.Close();

    std::cout << "Wrote with " << written << " DIGI-only events in total.\n";
    return written;
}

void printUsage(const char* argv0)
{
    std::cerr
        << "Usage:\n"
        << "  " << argv0 << " digi_path digi_run_id output_path [recreate_outputs]\n\n"
        << "Examples:\n"
        << "  " << argv0 << " /path/to/digi/splitted 1778243915 /path/to/output 1\n"
        << "  " << argv0 << " /path/to/digi/file.root 1778243915 /path/to/output 0\n\n"
        << "Output files are named XXXX_DIGI_ID_SPLIT.root, for example\n"
        << "XXXX_1778243915_0000000000.root.\n"
        << "recreate_outputs is optional: 1 deletes existing output files first "
        << "(default), 0 keeps and skips them.\n"
        << "If digi_path is a directory, only digi_run_id_*.root files are processed.\n";
}

} // namespace

int main(int argc, char** argv)
{
    if (argc != 4 && argc != 5) {
        printUsage(argv[0]);
        return 1;
    }

    try {
        const std::string digi_path = argv[1];
        const std::string digi_run_id = argv[2];
        const std::string output_path = argv[3];
        const bool recreate_outputs = parseRecreateOutputs(argc, argv);
        const auto digi_files = expandDigiInputs(digi_path, digi_run_id);

        if (digi_files.empty()) {
            throw std::runtime_error(
                "No digi ROOT files found for run " + digi_run_id + " in " + digi_path);
        }
        fs::create_directories(output_path);

        Long64_t total = 0;
        for (const auto& digi_file : digi_files) {
            total += buildOneFile(digi_file, digi_run_id, output_path, recreate_outputs);
        }

        std::cout << "All done. Total events written: " << total << '\n';
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << '\n';
        return 2;
    }

    return 0;
}
