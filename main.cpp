#include "wmFSPBWT.h"
#include <getopt.h>
#include <fstream>
#include <iostream>
#include <string>
#include <cmath>
#include <cstdint>

void printHelp(const char* programName) {
    std::cout << "Usage: " << programName << " [options]\n"
              << "  -i <file>     Input panel file (default: panel)\n"
              << "  -o <file>     Output file (default: <input>.out)\n"
              << "  -s <file>     Statistics info file (default: <output>.inf, customizable)\n"
              << "  -q <file>     External query sequence file (required for out-panel mode)\n"
              << "  -m <in/out>   Query mode: 'in' (in-panel, default) or 'out' (out-panel)\n"
              << "  -B <int>      Block size B (default: 64, supports 64 or 128)\n"
              << "  -f <int>      F value (default: 1, supports 1,2,3,4)\n"
              << "  -l <int>      Query length L (default: 1800, must be > 2*B-2)\n"
              << "  -h            Show this help message\n\n"
              << "Examples:\n"
              << "  " << programName << " -i panel.macs -o result.out -s result.inf -m in -l 2000\n"
              << "  " << programName << " -i data.macs -m out -q query.macs -B 128 -l 3000\n"
              << "  " << programName << " -i panel.macs -o myresult.out        # Automatically generates myresult.out and myresult.out.inf\n";
}

bool validateFiles(const std::string& inputFile, const std::string& outputFile,
                   const std::string& statsFile, const std::string& queryFile,
                   const std::string& queryMode) {
    // Input panel file
    {
        std::ifstream in(inputFile);
        if (!in.good()) {
            std::cerr << "Error: Input panel file '" << inputFile << "' does not exist or cannot be read.\n";
            return false;
        }
    }

    // Output file writability
    {
        std::ofstream out(outputFile, std::ios::app);
        if (!out.good()) {
            std::cerr << "Error: Cannot write to output file '" << outputFile << "'.\n";
            return false;
        }
    }

    // Stats file writability (only check if filename is determined)
    if (!statsFile.empty()) {
        std::ofstream sout(statsFile, std::ios::app);
        if (!sout.good()) {
            std::cerr << "Error: Cannot write to statistics file '" << statsFile << "'.\n";
            return false;
        }
    }

    // Query file for out-panel mode
    if (queryMode == "out") {
        if (queryFile.empty()) {
            std::cerr << "Error: Out-panel query mode (-m out) requires a query file (-q).\n";
            return false;
        }
        std::ifstream qf(queryFile);
        if (!qf.good()) {
            std::cerr << "Error: Query file '" << queryFile << "' does not exist or cannot be read.\n";
            return false;
        }
    }
    return true;
}

int main(int argc, char* argv[]) {
    std::string inputFile   = "panel";
    std::string outputFile  = "";
    std::string statsFile   = "";          // Empty means auto-generate
    std::string queryFile   = "";
    std::string queryMode   = "in";
    int B                   = 64;
    int F                   = 1;
    int queryLength         = 1800;

    int opt;
    const char* const short_opts = "i:I:o:O:s:S:q:Q:m:M:b:B:f:F:l:L:h:H";


    while ((opt = getopt(argc, argv, short_opts)) != -1) {
        switch (opt) {
        case 'i': case 'I':  // 大小写合并
            inputFile = optarg;
            break;

        case 'o': case 'O':
            outputFile = optarg;
            break;

        case 's': case 'S':
            statsFile = optarg;
            break;

        case 'q': case 'Q':
            queryFile = optarg;
            break;

        case 'm': case 'M':
            queryMode = optarg;
            if (queryMode != "in" && queryMode != "out") {
                std::cerr << "Error: Query mode must be 'in' or 'out'.\n";
                return 1;
            }
            break;

        case 'B': case 'b': // 重点：这里 B 和 b 都会进入同一个逻辑
            B = std::stoi(optarg);
            if (B != 64 && B != 128) {
                std::cerr << "Error: Currently only B = 64 or 128 is supported.\n";
                return 1;
            }
            break;

        case 'f': case 'F':
            F = std::stoi(optarg);
            if (F < 0) {
                std::cerr << "Error: F value must be >= 0.\n";
                return 1;
            }
            break;

        case 'l': case 'L':
            queryLength = std::stoi(optarg);
            break;

        case 'h': case 'H':
            printHelp(argv[0]);
            return 0;

        default:
            std::cerr << "Unknown option, use -h for help.\n";
            return 1;
        }
    }

    // Auto-generate output filename
    if (outputFile.empty()) {
        outputFile = inputFile + ".out";
    }

    // Auto-generate stats filename (if not specified by -s)
    if (statsFile.empty()) {
        statsFile = outputFile + ".inf";
    }

    // Parameter validity check
    // Constraint: L must be strictly greater than 2*B-2 (effectively >= 2*B-1)
    if (queryLength < 2 * B - 1) {
        std::cerr << "Error: Query length L must be strictly greater than 2*B-2 (Current 2*B = " << 2 * B << ", L = " << queryLength << ").\n";
        return 1;
    }

    if (!validateFiles(inputFile, outputFile, statsFile, queryFile, queryMode)) {
        return 1;
    }

    int ret = 0;

    if (B == 64) {
        wmFSPBWT<unsigned long long> CRY;
        CRY.F = F;
        CRY.B = B;
        CRY.T = std::pow(2, F);
        CRY.minSiteL = B * 2 - 1;

        ret = CRY.readMacsPanel(inputFile);
        std::cout << "Panel reading completed: " << ret << "\n";
        if (ret != 0) return ret;

        ret = CRY.makePanel();
        std::cout << "Fuzzy panel generation: " << ret << "\n";
        if (ret != 0) return ret;

        if (queryMode == "in") {
            ret = CRY.inPanelLongMatchQuery(queryLength, outputFile);
            std::cout << "In-panel query completed: " << ret << "\n";
        } else {
            ret = CRY.readMacsQuery(queryFile);
            std::cout << "Query sequence reading completed: " << ret << "\n";
            if (ret != 0) return ret;
            ret = CRY.outPanelLongMatchQuery(queryLength, outputFile);
            std::cout << "Out-panel query completed: " << ret << "\n";
        }

        CRY.outputInformationToFile(statsFile, queryMode);

    } else if (B == 128) {
        wmFSPBWT<__uint128_t> CRY;
        CRY.F = F;
        CRY.B = B;
        CRY.T = std::pow(2, F);
        CRY.minSiteL = B * 2 - 1;

        ret = CRY.readMacsPanel(inputFile);
        std::cout << "Panel reading completed: " << ret << "\n";
        if (ret != 0) return ret;

        ret = CRY.makePanel();
        std::cout << "Fuzzy panel generation: " << ret << "\n";
        if (ret != 0) return ret;

        if (queryMode == "in") {
            ret = CRY.inPanelLongMatchQuery(queryLength, outputFile);
            std::cout << "In-panel query completed: " << ret << "\n";
        } else {
            ret = CRY.readMacsQuery(queryFile);
            std::cout << "Query sequence reading completed: " << ret << "\n";
            if (ret != 0) return ret;
            ret = CRY.outPanelLongMatchQuery(queryLength, outputFile);
            std::cout << "Out-panel query completed: " << ret << "\n";
        }

        CRY.outputInformationToFile(statsFile, queryMode);
    }

    std::cout << "Query results saved to: " << outputFile << "\n";
    std::cout << "Statistics info saved to: " << statsFile << "\n";
    return 0;
}