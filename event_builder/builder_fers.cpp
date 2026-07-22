// FERS-only event builder for runs without a digitiser counterpart in the FCC TileCal beamtest data
// by M. Soldani, 2026 - developed with OpenAI Codex (GPT-5.5)

/* compilation:
g++ -std=c++17 builder_fers.cpp $(root-config --cflags --libs) -o builder_fers
*/

// usage: ./builder_fers <fers_root_file> <output_dir> [recreate_outputs=1]
/* --> examples:
./builder_fers  /home/msoldani/25-27_cern/26_05_bt_t10/data_prep/data/TEST_root/fers_root/merged/test_align_TStamp_0/Run1141.dat.root  ../data/TEST_root/global_root/splitted  1
./builder_fers  /eos/experiment/newtile/beamtests/26_05_t10/fers_root/merged/test_align_TStamp_0/Run1141.dat.root   /eos/experiment/newtile/beamtests/26_05_t10/global_root/splitted  1
*/

#include <TDirectory.h>
#include <TError.h>
#include <TFile.h>
#include <TKey.h>
#include <TObject.h>
#include <TString.h>
#include <TTimeStamp.h>
#include <TTree.h>

#include <algorithm>
#include <cctype>
#include <cstdlib>
#include <filesystem>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <streambuf>
#include <string>
#include <system_error>
#include <utility>
#include <vector>

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

TTree* treeOrNull(TFile& file, const std::string& tree_name)
{
    return dynamic_cast<TTree*>(file.Get(tree_name.c_str()));
}

TTree* requireDataTree(TFile& file)
{
    if (auto* tree = treeOrNull(file, "datas")) return tree;
    return firstTree(file);
}

void setFersBranchAddresses(TTree& fers_tree, FersBuffers& buffers)
{
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

void writeInfoTree(TTree* input_info)
{
    InfoBuffers buffers;

    if (input_info && input_info->GetEntries() > 0) {
        readInfoEntry(*input_info, buffers);
    }
    buffers.digi_run = 0;

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

bool hasFersBoard(const FersBuffers& buffers, int board)
{
    return buffers.Trg_Id[board] != 0;
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

std::string outputFilePath(const std::string& output_path, const std::string& fers_id)
{
    return (fs::path(output_path) / (fers_id + "_XXXXXXXXXX_0000000000.root")).string();
}

bool parseRecreateOutputs(int argc, char** argv)
{
    if (argc == 3) return true;

    const std::string value = argv[3];
    if (value == "0") return false;
    if (value == "1") return true;
    throw std::runtime_error("recreate_outputs must be 0 or 1");
}

Long64_t buildFersOnlyFile(
    const std::string& fers_path,
    const std::string& output_path,
    bool recreate_outputs)
{
    const std::string fers_id = firstDigitGroup(fs::path(fers_path).filename().string());
    if (fers_id.empty()) {
        throw std::runtime_error("Cannot extract FERS ID from ROOT filename: " + fers_path);
    }

    fs::create_directories(output_path);
    const std::string out_path = outputFilePath(output_path, fers_id);
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

    TFile fers_file(fers_path.c_str(), "READ");
    if (fers_file.IsZombie()) throw std::runtime_error("Cannot open FERS ROOT file: " + fers_path);
    TTree* fers_tree = requireDataTree(fers_file);
    TTree* info_tree = treeOrNull(fers_file, "info");
    if (!info_tree) {
        std::cerr << "Warning: no FERS info tree found in " << fers_path << ".\n";
    }

    TFile out_file(out_path.c_str(), "CREATE");
    if (out_file.IsZombie()) throw std::runtime_error("Cannot create output ROOT file: " + out_path);

    DigiBuffers digi_output;
    FersBuffers fers_input;
    FersBuffers fers_output;
    setFersBranchAddresses(*fers_tree, fers_input);
    auto* out_digi = createDigiOutputTree(digi_output);
    auto* out_fers = createFersOutputTree(fers_output);

    const Long64_t entries = fers_tree->GetEntries();
    Long64_t written = 0;
    Long64_t skipped_empty = 0;
    for (Long64_t entry = 0; entry < entries; ++entry) {
        const Long64_t processed = entry + 1;

        fers_tree->GetEntry(entry);
        int active_boards = 0;
        for (int board = 0; board < NFERSBDS; ++board) {
            if (hasFersBoard(fers_input, board)) ++active_boards;
        }
        if (active_boards == 0) {
            ++skipped_empty;
            if (processed % 10000 == 0 || processed == entries) {
                std::cout << "Processed " << processed << "/" << entries
                          << " input FERS entries; in total: wrote " << written
                          << ", skipped empty " << skipped_empty << ".\n";
            }
            continue;
        }

        for (int board = 0; board < NFERSBDS; ++board) {
            if (!hasFersBoard(fers_input, board)) continue;

            resetDigi(digi_output);
            resetFers(fers_output);
            copyFersBoard(fers_input, fers_output, board);
            digi_output.is_fers[board] = 1;

            out_digi->Fill();
            out_fers->Fill();
            ++written;
        }

        if (processed % 10000 == 0 || processed == entries) {
            std::cout << "Processed " << processed << "/" << entries
                      << " input FERS entries; in total: wrote " << written
                      << ", skipped empty " << skipped_empty << ".\n";
        }
    }

    out_file.cd();
    out_digi->Write();
    out_fers->Write();
    writeInfoTree(info_tree);
    out_file.Close();

    std::cout << "All done: written with " << written << " FERS-only events in total";
    if (skipped_empty > 0) {
        std::cout << " (" << skipped_empty << " empty input entries skipped)";
    }
    std::cout << ".\n";
    return written;
}

void printUsage(const char* argv0)
{
    std::cerr
        << "Usage:\n"
        << "  " << argv0 << " fers.root output_path [recreate_outputs]\n\n"
        << "Examples:\n"
        << "  " << argv0 << " Run1141.dat.root /path/to/output 1\n\n"
        << "Output files are named FERS_ID_XXXXXXXXXX_0000000000.root, for example\n"
        << "1141_XXXXXXXXXX_0000000000.root.\n"
        << "recreate_outputs is optional: 1 deletes existing output files first "
        << "(default), 0 keeps and skips them.\n";
}

} // namespace

int main(int argc, char** argv)
{
    if (argc != 3 && argc != 4) {
        printUsage(argv[0]);
        return 1;
    }

    std::string log_prefix;
    try {
        const std::string fers_path = argv[1];
        const std::string output_path = argv[2];
        const bool recreate_outputs = parseRecreateOutputs(argc, argv);
        const std::string fers_id = firstDigitGroup(fs::path(fers_path).filename().string());
        log_prefix = "[" + (fers_id.empty() ? "?" : fers_id) + "] ";
        ScopedStreamPrefix cout_prefix(std::cout, log_prefix);
        ScopedStreamPrefix cerr_prefix(std::cerr, log_prefix);
        SetErrorHandler(rootErrorHandler);
        buildFersOnlyFile(fers_path, output_path, recreate_outputs);
    } catch (const std::exception& e) {
        std::cerr << log_prefix << "Error: " << e.what() << '\n';
        return 2;
    }

    return 0;
}
