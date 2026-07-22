// event builder for the FCC TileCal beamtest data
// by M. Soldani, 2026 - developed with OpenAI Codex (GPT-5.5)

/* compilation:
g++ -std=c++17 builder.cpp $(root-config --cflags --libs) -o builder
*/

// usage: ./builder <sync_list.txt> <fers_root_file> <digi_root_dir> <run_id> <output_dir> [recreate_outputs=1]
/* --> examples:
./builder  ../data/output/sync_lists/1059.txt  /home/msoldani/25-27_cern/26_05_bt_t10/data_prep/data/TEST_root/fers_root/merged/test_align_TStamp_0/Run1059.dat.root  /home/msoldani/25-27_cern/26_05_bt_t10/data_prep/data/TEST_root/digi_root/splitted  1778258908  ../data/TEST_root/global_root/splitted  1
./builder  ../data/output/sync_lists/1059.txt  /eos/experiment/newtile/beamtests/26_05_t10/fers_root/merged/test_align_TStamp_0/Run1059.dat.root  /eos/experiment/newtile/beamtests/26_05_t10/digi_root/splitted  1778258908  /eos/experiment/newtile/beamtests/26_05_t10/global_root/splitted  1
*/

#include <TDirectory.h>
#include <TError.h>
#include <TFile.h>
#include <TKey.h>
#include <TLeaf.h>
#include <TObject.h>
#include <TString.h>
#include <TTimeStamp.h>
#include <TTree.h>
#include <TTreeFormula.h>

#include <algorithm>
#include <array>
#include <cctype>
#include <cstdlib>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <map>
#include <memory>
#include <streambuf>
#include <sstream>
#include <stdexcept>
#include <string>
#include <system_error>
#include <unordered_set>
#include <unordered_map>
#include <utility>
#include <vector>

#include <sys/wait.h>
#include <unistd.h>

#define NWFSAMPLES 1030
#define NFERSBDS 2 
#define NFERSCHS1BD 64

namespace fs = std::filesystem;

namespace {

class LinePrefixBuf : public std::streambuf {
public:
    LinePrefixBuf(std::streambuf* dest, std::string prefix)
        : dest_(dest), prefix_(std::move(prefix))
    {
    }

protected:
    int overflow(int ch) override
    {
        if (ch == traits_type::eof()) return traits_type::not_eof(ch);
        if (at_line_start_) {
            dest_->sputn(prefix_.data(), static_cast<std::streamsize>(prefix_.size()));
            at_line_start_ = false;
        }
        if (dest_->sputc(static_cast<char>(ch)) == traits_type::eof()) return traits_type::eof();
        if (ch == '\n') at_line_start_ = true;
        return ch;
    }

    int sync() override
    {
        return dest_->pubsync();
    }

private:
    std::streambuf* dest_ = nullptr;
    std::string prefix_;
    bool at_line_start_ = true;
};

class ScopedStreamPrefix {
public:
    ScopedStreamPrefix(std::ostream& stream, const std::string& prefix)
        : stream_(stream), old_buf_(stream.rdbuf()), prefix_buf_(old_buf_, prefix)
    {
        stream_.rdbuf(&prefix_buf_);
    }

    ~ScopedStreamPrefix()
    {
        stream_.rdbuf(old_buf_);
    }

private:
    std::ostream& stream_;
    std::streambuf* old_buf_ = nullptr;
    LinePrefixBuf prefix_buf_;
};

void rootErrorHandler(int level, Bool_t abort_bool, const char* location, const char* msg)
{
    const char* label = "Info";
    if (level >= kFatal) {
        label = "Fatal";
    } else if (level >= kSysError) {
        label = "SysError";
    } else if (level >= kError) {
        label = "Error";
    } else if (level >= kWarning) {
        label = "Warning";
    }

    std::cerr << label << " in <" << location << ">: " << msg << '\n';
    if (abort_bool) std::abort();
}

struct SyncRow {
    long long digi = -1;
    std::array<long long, NFERSBDS> fers;

    SyncRow()
    {
        fers.fill(-1);
    }
};

struct FersIndex {
    std::array<std::unordered_map<long long, Long64_t>, NFERSBDS> boards;
};

struct FersEntries {
    std::array<Long64_t, NFERSBDS> boards;

    FersEntries()
    {
        boards.fill(-1);
    }
};

struct OutputRecord {
    long long digi_id = -1;
    long long fers_id = -1;
    Long64_t digi_entry = -1;
    FersEntries fers_entries;
};

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

std::vector<long long> parseNumbers(std::string line)
{
    for (char& c : line) {
        if (c == ',' || c == '\t') c = ' ';
    }

    std::stringstream ss(line);
    std::vector<long long> values;
    long long value = 0;
    while (ss >> value) values.push_back(value);
    return values;
}

std::vector<SyncRow> readSyncRows(const std::string& path)
{
    std::ifstream in(path);
    if (!in) throw std::runtime_error("Cannot open sync list: " + path);

    std::vector<SyncRow> rows;
    std::string line;
    while (std::getline(in, line)) {
        if (line.empty() || line[0] == '#') continue;
        const auto values = parseNumbers(line);
        if (values.size() < static_cast<std::size_t>(1 + NFERSBDS)) continue;

        SyncRow row;
        row.digi = values[0];
        for (int board = 0; board < NFERSBDS; ++board) {
            row.fers[board] = values[1 + board];
        }
        rows.push_back(row);
    }

    if (rows.empty()) throw std::runtime_error("No usable rows in sync list: " + path);
    return rows;
}

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

TTree* treeOrNull(TFile& file, const std::string& tree_name)
{
    return dynamic_cast<TTree*>(file.Get(tree_name.c_str()));
}

TTree* requireDataTree(TFile& file)
{
    if (auto* tree = treeOrNull(file, "datas")) return tree;
    return firstTree(file);
}

void requireBranch(TTree& tree, const char* branch_name)
{
    if (!tree.GetBranch(branch_name)) {
        throw std::runtime_error(
            "Missing branch '" + std::string(branch_name) + "' in " + tree.GetName());
    }
}

TLeaf* requireLeaf(TTree& tree, const char* leaf_name)
{
    auto* leaf = tree.GetLeaf(leaf_name);
    if (!leaf) {
        throw std::runtime_error(
            "Missing leaf '" + std::string(leaf_name) + "' in " + tree.GetName());
    }
    return leaf;
}

void requireBranches(TTree& tree, const std::vector<const char*>& branch_names)
{
    for (const char* branch_name : branch_names) requireBranch(tree, branch_name);
}

void requireReadableEntry(TTree& tree, Long64_t entry)
{
    const Long64_t bytes = tree.GetEntry(entry);
    if (bytes <= 0) {
        throw std::runtime_error(
            "Cannot read entry " + std::to_string(entry) + " from tree " + tree.GetName());
    }
}

void checkAllEntriesReadable(TTree& tree)
{
    const Long64_t entries = tree.GetEntries();
    if (entries <= 0) {
        throw std::runtime_error("No entries found in tree " + std::string(tree.GetName()));
    }
    for (Long64_t entry = 0; entry < entries; ++entry) {
        requireReadableEntry(tree, entry);
    }
}

long long formulaValue(TTreeFormula& formula)
{
    formula.UpdateFormulaLeaves();
    return std::llround(formula.EvalInstance());
}

void insertIndex(
    std::unordered_map<long long, Long64_t>& index,
    std::unordered_set<long long>& duplicates,
    long long value,
    Long64_t entry,
    const std::string& label)
{
    if (value < 0) return;
    if (duplicates.find(value) != duplicates.end()) {
        std::cerr << "Warning: duplicate FERS " << label << " Trg_Id " << value
                  << "; ignoring entry " << entry << ".\n";
        return;
    }

    const auto [it, inserted] = index.emplace(value, entry);
    if (!inserted) {
        std::cerr << "Warning: duplicate FERS " << label << " Trg_Id " << value
                  << "; dropping entries " << it->second << " and " << entry << ".\n";
        index.erase(it);
        duplicates.insert(value);
    }
}

FersIndex buildFersIndex(TTree& fers_tree)
{
    FersIndex index;
    std::array<std::unordered_set<long long>, NFERSBDS> duplicates;
    std::vector<std::unique_ptr<TTreeFormula>> trg_formulas;
    trg_formulas.reserve(NFERSBDS);
    for (int board = 0; board < NFERSBDS; ++board) {
        trg_formulas.push_back(std::make_unique<TTreeFormula>(
            ("trg" + std::to_string(board)).c_str(),
            ("Trg_Id[" + std::to_string(board) + "]").c_str(),
            &fers_tree));
    }

    const Long64_t entries = fers_tree.GetEntries();
    for (Long64_t entry = 0; entry < entries; ++entry) {
        fers_tree.GetEntry(entry);
        std::array<long long, NFERSBDS> values = {};
        int active_board = -1;
        int active_boards = 0;
        for (int board = 0; board < NFERSBDS; ++board) {
            values[board] = formulaValue(*trg_formulas[board]);
            if (values[board] != 0) {
                active_board = board;
                ++active_boards;
            }
        }

        if (active_boards == 1) {
            insertIndex(
                index.boards[active_board],
                duplicates[active_board],
                values[active_board],
                entry,
                "board" + std::to_string(active_board));
        }
    }

    return index;
}

void setFersBranchAddresses(TTree& fers_tree, FersBuffers& buffers)
{
    requireBranches(fers_tree, {
        "TStamp_us", "Num_Hits", "dTRef", "Trg_Id", "ch_mask",
        "data_type", "PHA_LG", "PHA_HG", "ToA", "ToT"});

    fers_tree.SetBranchAddress("TStamp_us", buffers.TStamp_us);
    fers_tree.SetBranchAddress("Num_Hits", buffers.Num_Hits);
    fers_tree.SetBranchAddress("dTRef", buffers.dTRef);
    fers_tree.SetBranchAddress("Trg_Id", buffers.Trg_Id);
    fers_tree.SetBranchAddress("ch_mask", buffers.ch_mask);
    fers_tree.SetBranchAddress("data_type", buffers.data_type);
    fers_tree.SetBranchAddress("PHA_LG", buffers.PHA_LG);
    fers_tree.SetBranchAddress("PHA_HG", buffers.PHA_HG);
    fers_tree.SetBranchAddress("ToA", buffers.ToA);
    fers_tree.SetBranchAddress("ToT", buffers.ToT);
}

void setDigiBranchAddresses(TTree& digi_tree, DigiBuffers& buffers)
{
    requireBranches(digi_tree, {
        "run", "event0", "trigger_ts", "event_id", "ro_time",
        "wave0", "wave1", "wave2", "wave3", "wave4", "wave5", "wave6", "wave7"});

    const std::string trigger_type = requireLeaf(digi_tree, "trigger_ts")->GetTypeName();
    const std::string event_id_type = requireLeaf(digi_tree, "event_id")->GetTypeName();
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

long long digiEventId(const DigiBuffers& buffers)
{
    return buffers.event_id_is_long
        ? buffers.event_id_long
        : std::llround(buffers.event_id);
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

void readInfoEntry(TTree& info_tree, InfoBuffers& buffers)
{
    TString* file_format = nullptr;
    TString* janus_rel = nullptr;
    TString* acq_mode = nullptr;
    TTimeStamp* time_UTC = nullptr;
    TString* time_unit = nullptr;

    auto fersInfoName = [&](const char* name) {
        const std::string prefixed = std::string("FERS_") + name;
        return info_tree.GetBranch(prefixed.c_str()) ? prefixed : std::string(name);
    };

    info_tree.SetBranchAddress(fersInfoName("n_boards").c_str(), &buffers.n_boards);
    info_tree.SetBranchAddress(fersInfoName("n_channels").c_str(), &buffers.n_channels);
    info_tree.SetBranchAddress(fersInfoName("max_hits").c_str(), &buffers.max_hits);
    info_tree.SetBranchAddress(fersInfoName("board_mod").c_str(), &buffers.board_mod);
    info_tree.SetBranchAddress(fersInfoName("file_format").c_str(), &file_format);
    info_tree.SetBranchAddress(fersInfoName("janus_rel").c_str(), &janus_rel);
    info_tree.SetBranchAddress(fersInfoName("acq_mode").c_str(), &acq_mode);
    info_tree.SetBranchAddress(fersInfoName("run").c_str(), &buffers.run);
    info_tree.SetBranchAddress(fersInfoName("time_epoch").c_str(), &buffers.time_epoch);
    info_tree.SetBranchAddress(fersInfoName("time_UTC").c_str(), &time_UTC);
    info_tree.SetBranchAddress(fersInfoName("e_Nbins").c_str(), &buffers.e_Nbins);
    info_tree.SetBranchAddress(fersInfoName("time_LSB_ns").c_str(), &buffers.time_LSB_ns);
    info_tree.SetBranchAddress(fersInfoName("time_unit").c_str(), &time_unit);

    info_tree.GetEntry(0);

    if (file_format) buffers.file_format = *file_format;
    if (janus_rel) buffers.janus_rel = *janus_rel;
    if (acq_mode) buffers.acq_mode = *acq_mode;
    if (time_UTC) buffers.time_UTC = *time_UTC;
    if (time_unit) buffers.time_unit = *time_unit;
}

void writeInfoTree(TTree* input_info, ULong64_t digi_run)
{
    InfoBuffers buffers;
    buffers.digi_run = digi_run;

    if (input_info && input_info->GetEntries() > 0) {
        readInfoEntry(*input_info, buffers);
        buffers.digi_run = digi_run;
    }

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

void resetDigi(DigiBuffers& buffers)
{
    buffers.run = 0;
    buffers.event0 = -1;
    buffers.trigger_ts = 0;
    buffers.event_id = -1;
    buffers.trigger_ts_long = 0;
    buffers.event_id_long = -1;
    buffers.ro_time = 0;
    for (int i = 0; i < NWFSAMPLES; ++i) {
        buffers.wave0[i] = 0;
        buffers.wave1[i] = 0;
        buffers.wave2[i] = 0;
        buffers.wave3[i] = 0;
        buffers.wave4[i] = 0;
        buffers.wave5[i] = 0;
        buffers.wave6[i] = 0;
        buffers.wave7[i] = 0;
    }
    buffers.is_digi = 0;
    for (int board = 0; board < NFERSBDS; ++board) {
        buffers.is_fers[board] = 0;
    }
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

void copyFersBoard(const FersBuffers& input, FersBuffers& output, int board)
{
    output.TStamp_us[board] = input.TStamp_us[board];
    output.Num_Hits[board] = input.Num_Hits[board];
    output.dTRef[board] = input.dTRef[board];
    output.Trg_Id[board] = input.Trg_Id[board];
    output.ch_mask[board] = input.ch_mask[board];

    for (int ch = 0; ch < NFERSCHS1BD; ++ch) {
        output.data_type[board][ch] = input.data_type[board][ch];
        output.PHA_LG[board][ch] = input.PHA_LG[board][ch];
        output.PHA_HG[board][ch] = input.PHA_HG[board][ch];
        output.ToA[board][ch] = input.ToA[board][ch];
        output.ToT[board][ch] = input.ToT[board][ch];
    }
    output.is_fers[board] = 1;
}

std::map<long long, std::vector<const SyncRow*>> rowsByDigi(const std::vector<SyncRow>& rows)
{
    std::map<long long, std::vector<const SyncRow*>> out;
    for (const auto& row : rows) {
        if (row.digi < 0) continue;
        out[row.digi].push_back(&row);
    }
    return out;
}

std::unordered_set<long long> syncedFersIds(const std::vector<SyncRow>& rows)
{
    std::unordered_set<long long> out;
    for (const auto& row : rows) {
        for (int board = 0; board < NFERSBDS; ++board) {
            if (row.fers[board] >= 0) out.insert(row.fers[board]);
        }
    }
    return out;
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

std::string firstDigitGroup(const std::string& text)
{
    std::string digits;
    for (const char c : text) {
        if (std::isdigit(static_cast<unsigned char>(c))) {
            digits.push_back(c);
        } else if (!digits.empty()) {
            break;
        }
    }
    return digits;
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
    const std::string& fers_id,
    const std::string& digi_run_id,
    const std::string& digi_file_path)
{
    const std::string name = fers_id + "_" + digi_run_id + "_"
        + digiSplitId(digi_file_path, digi_run_id) + ".root";
    return (fs::path(output_path) / name).string();
}

bool parseRecreateOutputs(int argc, char** argv)
{
    if (argc == 6) return true;

    const std::string value = argv[6];
    if (value == "0") return false;
    if (value == "1") return true;
    throw std::runtime_error("recreate_outputs must be 0 or 1");
}

std::vector<std::string> filesToBuild(
    const std::vector<std::string>& digi_files,
    const std::string& output_path,
    const std::string& fers_id,
    const std::string& digi_run_id,
    bool recreate_outputs)
{
    if (recreate_outputs) return digi_files;

    std::vector<std::string> pending;
    for (const auto& digi_file : digi_files) {
        const std::string out_path = outputFilePath(output_path, fers_id, digi_run_id, digi_file);
        if (!fs::exists(out_path)) pending.push_back(digi_file);
    }
    return pending;
}

void checkFersRootFileContents(const std::string& fers_path)
{
    TFile fers_file(fers_path.c_str(), "READ");
    if (fers_file.IsZombie()) throw std::runtime_error("Cannot open FERS ROOT file: " + fers_path);
    TTree* fers_tree = requireDataTree(fers_file);
    FersBuffers buffers;
    setFersBranchAddresses(*fers_tree, buffers);
    checkAllEntriesReadable(*fers_tree);
}

void checkDigiRootFileContents(const std::string& digi_path)
{
    TFile digi_file(digi_path.c_str(), "READ");
    if (digi_file.IsZombie()) throw std::runtime_error("Cannot open digi ROOT file: " + digi_path);
    TTree* digi_tree = firstTree(digi_file);
    DigiBuffers buffers;
    setDigiBranchAddresses(*digi_tree, buffers);
    checkAllEntriesReadable(*digi_tree);
}

void checkRootFileInChild(const std::string& path, const std::string& label, void (*check)(const std::string&))
{
    const pid_t pid = fork();
    if (pid < 0) {
        throw std::runtime_error("Cannot fork while checking " + label + " ROOT file: " + path);
    }

    if (pid == 0) {
        try {
            check(path);
            _exit(0);
        } catch (const std::exception& e) {
            std::cerr << "Error while checking " << label << " ROOT file " << path
                      << ": " << e.what() << '\n';
            std::cerr.flush();
            _exit(2);
        } catch (...) {
            std::cerr << "Unknown error while checking " << label << " ROOT file "
                      << path << ".\n";
            std::cerr.flush();
            _exit(2);
        }
    }

    int status = 0;
    if (waitpid(pid, &status, 0) < 0) {
        throw std::runtime_error("Cannot wait for " + label + " ROOT file check: " + path);
    }
    if (WIFSIGNALED(status)) {
        throw std::runtime_error(
            label + " ROOT file check crashed with signal "
            + std::to_string(WTERMSIG(status)) + ": " + path);
    }
    if (!WIFEXITED(status) || WEXITSTATUS(status) != 0) {
        throw std::runtime_error(label + " ROOT file check failed: " + path);
    }
}

void checkFersRootFile(const std::string& fers_path)
{
    checkRootFileInChild(fers_path, "FERS", checkFersRootFileContents);
}

void checkDigiRootFiles(const std::vector<std::string>& digi_files)
{
    for (const auto& digi_path : digi_files) {
        checkRootFileInChild(digi_path, "digi", checkDigiRootFileContents);
    }
}

Long64_t buildOneFile(
    const std::string& digi_path,
    const std::string& fers_id,
    const std::string& digi_run_id,
    TTree& fers_tree,
    TTree* info_tree,
    const FersIndex& fers_index,
    const std::map<long long, std::vector<const SyncRow*>>& sync_by_digi,
    const std::unordered_set<long long>& sync_fers_ids,
    const std::string& output_path,
    bool recreate_outputs)
{
    const std::string out_path = outputFilePath(output_path, fers_id, digi_run_id, digi_path);
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
    setDigiBranchAddresses(*digi_tree, digi_input);
    auto* out_digi = createDigiOutputTree(digi_output);

    FersBuffers fers_input;
    FersBuffers fers_output;
    setFersBranchAddresses(fers_tree, fers_input);
    auto* out_fers = createFersOutputTree(fers_output);

    std::map<long long, Long64_t> digi_index;
    Long64_t min_digi_id = 0;
    Long64_t max_digi_id = -1;
    const Long64_t digi_entries = digi_tree->GetEntries();
    for (Long64_t entry = 0; entry < digi_entries; ++entry) {
        digi_tree->GetEntry(entry);
        const long long event_id = digiEventId(digi_input);
        digi_index.emplace(event_id, entry);
        if (max_digi_id < min_digi_id) {
            min_digi_id = event_id;
            max_digi_id = event_id;
        } else {
            min_digi_id = std::min<Long64_t>(min_digi_id, event_id);
            max_digi_id = std::max<Long64_t>(max_digi_id, event_id);
        }
    }

    std::vector<OutputRecord> records;
    Long64_t min_fers_entry = -1;
    Long64_t max_fers_entry = -1;

    for (const auto& [digi_value, digi_entry] : digi_index) {
        const auto rows_it = sync_by_digi.find(digi_value);
        if (rows_it == sync_by_digi.end()) {
            OutputRecord record;
            record.digi_id = digi_value;
            record.digi_entry = digi_entry;
            records.push_back(record);
            continue;
        }

        for (const SyncRow* row : rows_it->second) {
            OutputRecord record;
            record.digi_id = digi_value;
            record.digi_entry = digi_entry;
            bool missing_requested_fers = false;

            for (int board = 0; board < NFERSBDS; ++board) {
                const long long fers_id = row->fers[board];
                if (fers_id < 0) continue;

                if (record.fers_id < 0) record.fers_id = fers_id;
                const auto fers_it = fers_index.boards[board].find(fers_id);
                if (fers_it != fers_index.boards[board].end()) {
                    record.fers_entries.boards[board] = fers_it->second;
                    min_fers_entry = min_fers_entry < 0
                        ? fers_it->second
                        : std::min(min_fers_entry, fers_it->second);
                    max_fers_entry = std::max(max_fers_entry, fers_it->second);
                    continue;
                }

                std::cerr << "Warning: missing FERS board" << board << " match for digi "
                          << row->digi << " (fers" << board << " " << fers_id
                          << "). Dropping the row.\n";
                missing_requested_fers = true;
            }

            if (missing_requested_fers) continue;
            records.push_back(record);
        }
    }

    std::map<std::pair<long long, int>, Long64_t> extra_fers;
    for (int board = 0; board < NFERSBDS; ++board) {
        for (const auto& [id, entry] : fers_index.boards[board]) {
            if (sync_fers_ids.find(id) != sync_fers_ids.end()) continue;
            extra_fers[{id, board}] = entry;
        }
    }

    for (const auto& [key, entry] : extra_fers) {
        if (min_fers_entry < 0) continue;
        if (entry < min_fers_entry || entry > max_fers_entry) continue;

        OutputRecord record;
        const auto [id, board] = key;
        record.fers_id = id;
        record.digi_entry = -1;
        record.fers_entries.boards[board] = entry;
        records.push_back(record);
    }

    std::vector<long long> digi_anchors;
    std::vector<long long> fers_anchors;
    for (const auto& record : records) {
        if (record.digi_id < 0 || record.fers_id < 0) continue;
        digi_anchors.push_back(record.digi_id);
        fers_anchors.push_back(record.fers_id);
    }
    std::sort(digi_anchors.begin(), digi_anchors.end());
    std::sort(fers_anchors.begin(), fers_anchors.end());

    auto orderGroup = [&](const OutputRecord& record) {
        if (record.digi_id >= 0 && record.fers_id >= 0) {
            const auto it = std::lower_bound(digi_anchors.begin(), digi_anchors.end(), record.digi_id);
            return 2 * static_cast<long long>(it - digi_anchors.begin());
        }
        if (record.digi_id >= 0) {
            const auto it = std::lower_bound(digi_anchors.begin(), digi_anchors.end(), record.digi_id);
            return 2 * static_cast<long long>(it - digi_anchors.begin()) - 1;
        }
        const auto it = std::lower_bound(fers_anchors.begin(), fers_anchors.end(), record.fers_id);
        return 2 * static_cast<long long>(it - fers_anchors.begin()) - 1;
    };

    std::stable_sort(records.begin(), records.end(), [&](const OutputRecord& a, const OutputRecord& b) {
        const long long group_a = orderGroup(a);
        const long long group_b = orderGroup(b);
        if (group_a != group_b) return group_a < group_b;

        if (a.digi_id >= 0 && b.digi_id >= 0 && a.digi_id != b.digi_id) {
            return a.digi_id < b.digi_id;
        }
        if (a.fers_id >= 0 && b.fers_id >= 0 && a.fers_id != b.fers_id) {
            return a.fers_id < b.fers_id;
        }
        if ((a.digi_id >= 0) != (b.digi_id >= 0)) return a.digi_id >= 0;
        return false;
    });

    Long64_t written = 0;
    for (const auto& record : records) {
        if (record.digi_entry >= 0) {
            digi_tree->GetEntry(record.digi_entry);
            copyDigi(digi_input, digi_output);
        } else {
            resetDigi(digi_output);
        }

        resetFers(fers_output);
        for (int board = 0; board < NFERSBDS; ++board) {
            if (record.fers_entries.boards[board] < 0) continue;
            fers_tree.GetEntry(record.fers_entries.boards[board]);
            copyFersBoard(fers_input, fers_output, board);
        }

        for (int board = 0; board < NFERSBDS; ++board) {
            digi_output.is_fers[board] = fers_output.is_fers[board];
        }
        fers_output.is_digi = digi_output.is_digi;

        out_digi->Fill();
        out_fers->Fill();
        ++written;
    }

    out_file.cd();
    out_digi->Write();
    out_fers->Write();

    out_file.cd();
    writeInfoTree(info_tree, std::stoull(digi_run_id));

    out_file.Close();

    std::cout << "Wrote with " << written << " events in total.\n";
    return written;
}

void printUsage(const char* argv0)
{
    std::cerr
        << "Usage:\n"
        << "  " << argv0 << " sync_list.txt fers.root digi_path digi_run_id output_path [recreate_outputs]\n\n"
        << "Examples:\n"
        << "  " << argv0 << " sync.txt Run1059.dat.root /path/to/digi/splitted 1778258908 /path/to/output 1\n"
        << "  " << argv0 << " sync.txt Run1059.dat.root /path/to/digi/file.root 1778258908 /path/to/output 0\n\n"
        << "Output files are named FERS_ID_DIGI_ID_SPLIT.root, for example\n"
        << "1059_1778258908_0000000000.root.\n"
        << "recreate_outputs is optional: 1 deletes existing output files first "
        << "(default), 0 keeps and skips them.\n"
        << "The sync list must have columns: digi_event_id followed by one FERS Trg_Id "
        << "column per board.\n"
        << "If digi_path is a directory, only digi_run_id_*.root files are processed.\n"
        << "The digi tree is matched with branch/expression event_id.\n"
        << "The FERS tree is indexed by the single non-zero Trg_Id board in each entry.\n";
}

} // namespace

int main(int argc, char** argv)
{
    if (argc != 6 && argc != 7) {
        printUsage(argv[0]);
        return 1;
    }

    std::string log_prefix;
    try {
        const std::string sync_path = argv[1];
        const std::string fers_path = argv[2];
        const std::string digi_path = argv[3];
        const std::string digi_run_id = argv[4];
        const std::string output_path = argv[5];
        const bool recreate_outputs = parseRecreateOutputs(argc, argv);
        const std::string fers_id = firstDigitGroup(fs::path(fers_path).filename().string());
        if (fers_id.empty()) {
            throw std::runtime_error("Cannot extract FERS ID from ROOT filename: " + fers_path);
        }
        log_prefix = "[" + fers_id + ", " + digi_run_id + ", "
            + digiSplitId(digi_path, digi_run_id) + "] ";
        ScopedStreamPrefix cout_prefix(std::cout, log_prefix);
        ScopedStreamPrefix cerr_prefix(std::cerr, log_prefix);
        SetErrorHandler(rootErrorHandler);

        const auto digi_files = expandDigiInputs(digi_path, digi_run_id);

        if (digi_files.empty()) {
            throw std::runtime_error(
                "No digi ROOT files found for run " + digi_run_id + " in " + digi_path);
        }
        fs::create_directories(output_path);
        const auto pending_digi_files = filesToBuild(
            digi_files, output_path, fers_id, digi_run_id, recreate_outputs);

        if (pending_digi_files.empty()) {
            for (const auto& digi_file : digi_files) {
                const std::string out_path = outputFilePath(output_path, fers_id, digi_run_id, digi_file);
                std::cout << "Keeping existing " << out_path << "; skipping.\n";
            }
            std::cout << "All done. Total events written: 0\n";
            return 0;
        }

        checkFersRootFile(fers_path);
        checkDigiRootFiles(pending_digi_files);

        const auto sync_rows = readSyncRows(sync_path);
        const auto sync_by_digi = rowsByDigi(sync_rows);
        const auto sync_fers_ids = syncedFersIds(sync_rows);

        TFile fers_file(fers_path.c_str(), "READ");
        if (fers_file.IsZombie()) throw std::runtime_error("Cannot open FERS ROOT file: " + fers_path);
        TTree* fers_tree = requireDataTree(fers_file);
        TTree* info_tree = treeOrNull(fers_file, "info");
        if (!info_tree) {
            std::cerr << "Warning: no FERS info tree found in " << fers_path << ".\n";
        }
        const FersIndex fers_index = buildFersIndex(*fers_tree);

        std::cout << "Indexed FERS entries:";
        for (int board = 0; board < NFERSBDS; ++board) {
            std::cout << " board" << board << "=" << fers_index.boards[board].size();
        }
        std::cout << '\n';

        Long64_t total = 0;
        for (const auto& digi_path : digi_files) {
            total += buildOneFile(
                digi_path, fers_id, digi_run_id, *fers_tree, info_tree,
                fers_index, sync_by_digi, sync_fers_ids, output_path, recreate_outputs);
        }

        std::cout << "All done. Total events written: " << total << '\n';
    } catch (const std::exception& e) {
        std::cerr << log_prefix << "Error: " << e.what() << '\n';
        return 2;
    }

    return 0;
}
