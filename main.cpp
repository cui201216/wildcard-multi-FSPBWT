#include "wmFSPBWT.h"
#include <getopt.h>
#include <fstream>
#include <iostream>
#include <string>
#include <cmath>
#include <cstdint>
#include <unistd.h>

using namespace std;

void printHelp(const char* programName) {
    cout << "Usage: " << programName << " [command] [options]\n\n";
    cout << "Commands:\n";
    cout << "  (no command)          Traditional mode: read panel, build index in memory, query directly\n";
    cout << "  build                 Build index and save to file\n";
    cout << "  self                  Load index and perform in-panel (self) query\n";
    cout << "  query                 Load index and perform out-panel query\n\n";
    cout << "Options (common):\n";
    cout << "  -i <file>             Input panel file (traditional/build) or index file (self/query)\n";
    cout << "  -o <file>             Output file (query result) or index file (build)\n";
    cout << "  -q <file>             Query file (only for out-panel / query mode)\n";
    cout << "  -m <in|out>           Query mode: 'in' or 'out' (traditional mode only)\n";
    cout << "  -l <int>              Minimum match length L (default: 1800)\n";
    cout << "  -B <64|128>           Block size B (default: 64)\n";
    cout << "  -f <1-4>              F value (default: 1)\n";
    cout << "  --macs                Input file is in MACS format\n";
    cout << "  -h                    Show this help\n\n";
    cout << "Examples:\n";
    cout << "  " << programName << " -i panel.vcf -m in -l 2000\n";
    cout << "  " << programName << " build -i panel.vcf -o panel.idx --macs -B 128 -f 2\n";
    cout << "  " << programName << " self -i panel.idx -o result.in.txt -l 1800\n";
    cout << "  " << programName << " query -i panel.idx -q query.vcf -o result.out.txt -l 1800\n";
}

int main(int argc, char* argv[]) {
    if (argc < 2) {
        printHelp(argv[0]);
        return 1;
    }

    string command;
    string panelFile, indexFile, queryFile, outputFile;
    string queryMode = "in";
    int B = 64;
    int F = 1;
    int L = 1800;
    bool isMacs = false;

    // 如果第一个参数不是选项（不是以 - 开头），则认为是子命令
    if (argv[1][0] != '-') {
        command = argv[1];
        optind = 2;  // 从第三个参数开始解析选项
    } else {
        command = "";  // 传统模式
        optind = 1;
    }

    // 统一选项解析
    int opt;
    const char* short_opts = "i:o:q:m:l:B:f:h";
    static struct option long_options[] = {
        {"macs", no_argument, 0, 'x'},
        {0, 0, 0, 0}
    };
    int option_index = 0;

    while ((opt = getopt_long(argc, argv, short_opts, long_options, &option_index)) != -1) {
        switch (opt) {
            case 'i': panelFile = optarg; break;
            case 'o': outputFile = optarg; break;
            case 'q': queryFile = optarg; break;
            case 'm':
                queryMode = optarg;
                if (queryMode != "in" && queryMode != "out") {
                    cerr << "Error: -m must be 'in' or 'out'\n";
                    return 1;
                }
                break;
            case 'l': L = stoi(optarg); break;
            case 'B': B = stoi(optarg); break;
            case 'f': F = stoi(optarg); break;
            case 'x': isMacs = true; break;
            case 'h': printHelp(argv[0]); return 0;
            default:
                cerr << "Unknown option\n";
                printHelp(argv[0]);
                return 1;
        }
    }

    // 参数校验
    if (L < 2 * B - 1) {
        cerr << "Error: L must be > 2*B-2\n";
        return 1;
    }

    int ret = 0;

    // 根据 B 选择 Syllable 类型
    if (B == 64) {
        wmFSPBWT<unsigned long long> CRY;
        CRY.B = B;
        CRY.F = F;
        CRY.T = (int)pow(2, F);
        CRY.minSiteL = B * 2 - 1;

        if (command == "build") {
            // build 模式：读面板 → 构建 → 保存索引
            if (panelFile.empty() || outputFile.empty()) {
                cerr << "build requires -i <panel> -o <index>\n";
                return 1;
            }
            cout << "[Build] Reading panel...\n";
            if (isMacs) ret = CRY.readMacsPanel(panelFile);
            else ret = CRY.readVcfPanel(panelFile);
            if (ret != 0) return ret;

            cout << "[Build] Making fuzzy panel...\n";
            ret = CRY.makePanel();
            if (ret != 0) return ret;

            cout << "[Build] Saving index to " << outputFile << "\n";
            CRY.save(outputFile);
            cout << "Index built and saved successfully.\n";
        }
        else if (command == "self" || command == "query") {
            // self / query 模式：加载索引 → 执行查询
            if (panelFile.empty() || outputFile.empty()) {
                cerr << command << " requires -i <index> -o <output>\n";
                return 1;
            }
            cout << "[Load] Loading index from " << panelFile << "\n";
            CRY.load(panelFile);

            if (command == "self") {
                cout << "[Query] Performing in-panel (self) query...\n";
                ret = CRY.inPanelLongMatchQuery(L, outputFile);
            } else {  // query
                if (queryFile.empty()) {
                    cerr << "query requires -q <query_file>\n";
                    return 1;
                }
                cout << "[Query] Reading query file...\n";
                if (isMacs) ret = CRY.readMacsQuery(queryFile);
                else ret = CRY.readVcfQuery(queryFile);
                if (ret != 0) return ret;

                cout << "[Query] Performing out-panel query...\n";
                ret = CRY.outPanelLongMatchQuery(L, outputFile);
            }
            if (ret == 0) {
                cout << "Query completed. Results saved to " << outputFile << "\n";
            }
        }
        else {
            // 传统模式（无 command）
            if (panelFile.empty()) {
                cerr << "No input panel specified (-i)\n";
                return 1;
            }
            cout << "[Traditional] Reading panel...\n";
            if (isMacs) ret = CRY.readMacsPanel(panelFile);
            else ret = CRY.readVcfPanel(panelFile);
            if (ret != 0) return ret;

            cout << "[Traditional] Making fuzzy panel...\n";
            ret = CRY.makePanel();
            if (ret != 0) return ret;

            if (queryMode == "in") {
                cout << "[Traditional] Performing in-panel query...\n";
                ret = CRY.inPanelLongMatchQuery(L, outputFile.empty() ? panelFile + ".out" : outputFile);
            } else {
                if (queryFile.empty()) {
                    cerr << "Out-panel mode requires -q <query_file>\n";
                    return 1;
                }
                cout << "[Traditional] Reading query...\n";
                if (isMacs) ret = CRY.readMacsQuery(queryFile);
                else ret = CRY.readVcfQuery(queryFile);
                if (ret != 0) return ret;

                cout << "[Traditional] Performing out-panel query...\n";
                ret = CRY.outPanelLongMatchQuery(L, outputFile.empty() ? panelFile + ".out" : outputFile);
            }
        }
    }
    else if (B == 128) {
        // 完全相同的逻辑，但使用 __uint128_t
        wmFSPBWT<__uint128_t> CRY;
        CRY.B = B;
        CRY.F = F;
        CRY.T = (int)pow(2, F);
        CRY.minSiteL = B * 2 - 1;

          if (command == "build") {
            // build 模式：读面板 → 构建 → 保存索引
            if (panelFile.empty() || outputFile.empty()) {
                cerr << "build requires -i <panel> -o <index>\n";
                return 1;
            }
            cout << "[Build] Reading panel...\n";
            if (isMacs) ret = CRY.readMacsPanel(panelFile);
            else ret = CRY.readVcfPanel(panelFile);
            if (ret != 0) return ret;

            cout << "[Build] Making fuzzy panel...\n";
            ret = CRY.makePanel();
            if (ret != 0) return ret;

            cout << "[Build] Saving index to " << outputFile << "\n";
            CRY.save(outputFile);
            cout << "Index built and saved successfully.\n";
        }
        else if (command == "self" || command == "query") {
            // self / query 模式：加载索引 → 执行查询
            if (panelFile.empty() || outputFile.empty()) {
                cerr << command << " requires -i <index> -o <output>\n";
                return 1;
            }
            cout << "[Load] Loading index from " << panelFile << "\n";
            CRY.load(panelFile);

            if (command == "self") {
                cout << "[Query] Performing in-panel (self) query...\n";
                ret = CRY.inPanelLongMatchQuery(L, outputFile);
            } else {  // query
                if (queryFile.empty()) {
                    cerr << "query requires -q <query_file>\n";
                    return 1;
                }
                cout << "[Query] Reading query file...\n";
                if (isMacs) ret = CRY.readMacsQuery(queryFile);
                else ret = CRY.readVcfQuery(queryFile);
                if (ret != 0) return ret;

                cout << "[Query] Performing out-panel query...\n";
                ret = CRY.outPanelLongMatchQuery(L, outputFile);
            }
            if (ret == 0) {
                cout << "Query completed. Results saved to " << outputFile << "\n";
            }
        }
        else {
            // 传统模式（无 command）
            if (panelFile.empty()) {
                cerr << "No input panel specified (-i)\n";
                return 1;
            }
            cout << "[Traditional] Reading panel...\n";
            if (isMacs) ret = CRY.readMacsPanel(panelFile);
            else ret = CRY.readVcfPanel(panelFile);
            if (ret != 0) return ret;

            cout << "[Traditional] Making fuzzy panel...\n";
            ret = CRY.makePanel();
            if (ret != 0) return ret;

            if (queryMode == "in") {
                cout << "[Traditional] Performing in-panel query...\n";
                ret = CRY.inPanelLongMatchQuery(L, outputFile.empty() ? panelFile + ".out" : outputFile);
            } else {
                if (queryFile.empty()) {
                    cerr << "Out-panel mode requires -q <query_file>\n";
                    return 1;
                }
                cout << "[Traditional] Reading query...\n";
                if (isMacs) ret = CRY.readMacsQuery(queryFile);
                else ret = CRY.readVcfQuery(queryFile);
                if (ret != 0) return ret;

                cout << "[Traditional] Performing out-panel query...\n";
                ret = CRY.outPanelLongMatchQuery(L, outputFile.empty() ? panelFile + ".out" : outputFile);
            }
        }
    }
    else {
        cerr << "Unsupported B value: " << B << "\n";
        return 1;
    }

    if (ret != 0) {
        cerr << "Operation failed with code: " << ret << endl;
    }

    return ret;
}