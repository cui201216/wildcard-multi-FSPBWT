#include "wmFSPBWT.h"

void printHelp(const char* programName) {
    std::cout << "用法: " << programName << " [-i input] [-o output] [-q panelFile] [-t in/out] [-B block_size] [-f F_value] [-l length]\n"
              << "  -i <file>  输入文件 (默认: out)\n"
              << "  -o <file>  输出文件 (默认: <input>.out)\n"
              << "  -q <file>  面板文件 (可选，用于面板外查询)\n"
              << "  -t <in/out> 查询模式：in（面板内查询，默认）或 out（面板外查询）\n"
              << "  -B <int>   音节大小 B (默认: 64)\n"
              << "  -f <int>   F 值 (默认: 1)\n"
              << "  -l <int>   查询长度 (默认: 1800)\n"
              << "示例: " << programName << " -i data.txt -q panel.txt -t out -B 128 -l 2000\n";
}

// 验证文件路径
bool validateFiles(const std::string& inputFile, const std::string& outputFile, const std::string& panelFile, const std::string& queryMode) {
    std::ifstream in(inputFile);
    if (!in.good()) {
        std::cerr << "错误: 输入文件 '" << inputFile << "' 不存在\n";
        return false;
    }
    in.close();
    std::ofstream out(outputFile, std::ios::app);
    if (!out.good()) {
        std::cerr << "错误: 无法写入输出文件 '" << outputFile << "'\n";
        return false;
    }
    out.close();
    if (queryMode == "out" && !panelFile.empty()) {
        std::ifstream panel(panelFile);
        if (!panel.good()) {
            std::cerr << "错误: 面板文件 '" << panelFile << "' 不存在\n";
            return false;
        }
        panel.close();
    } else if (queryMode == "out" && panelFile.empty()) {
        std::cerr << "错误: 面板外查询需要指定面板文件 (-q)\n";
        return false;
    }
    return true;
}

int main(int argc, char* argv[]) {
    std::string panelFile = "panel";    // 默认输入文件名
    std::string outputFile = "outTEST";           // 输出文件动态生成
    std::string queryFile= "query";            // 面板文件（可选）
    std::string queryMode = "in";     // 默认查询模式为面板内
    int B = 64;                       // 默认 B 值
    int F = 2;                        // 默认 F 值
    int queryLength = 200;           // 默认查询长度

    // 解析命令行参数
    int opt;
    while ((opt = getopt(argc, argv, "i:I:o:O:q:Q:t:T:b:B:f:F:l:L:h:H")) != -1) {
        try {
            switch (opt) {
            case 'i':
            case 'I':
                panelFile = optarg;
                break;
            case 'o':
            case 'O':
                outputFile = optarg;
                break;
            case 'q':
            case 'Q':
                queryFile = optarg;
                break;
            case 't':
            case 'T':
                queryMode = optarg;
                if (queryMode != "in" && queryMode != "out") {
                    std::cerr << "错误: 查询模式必须为 'in' 或 'out'\n";
                    return 1;
                }
                break;
            case 'b':
            case 'B':
                B = std::stoi(optarg);
                break;
            case 'f':
            case 'F':
                F = std::stoi(optarg);
                break;
            case 'l':
            case 'L':
                queryLength = std::stoi(optarg);
                break;
            case 'h':
            case 'H':
                printHelp(argv[0]);
                return 0;
            default:
                std::cerr << "错误: 未知选项，使用 -h 查看用法\n";
                return 1;
            }
        } catch (const std::exception& e) {
            std::cerr << "错误: 无效的参数值 '" << optarg << "'\n";
            return 1;
        }
    }

    // 自动生成输出文件名（若未指定）
    if (outputFile.empty()) {
        outputFile = panelFile + ".out";
    }

    // 验证参数
    if (B <= 0 || F < 0 || queryLength <= 0) {
        std::cerr << "错误: B、F 和查询长度必须为正数\n";
        return 1;
    }
    if (!validateFiles(panelFile, outputFile, queryFile, queryMode)) {
        return 1;
    }

    // // 输出参数信息
    // std::cout << "参数:\n"
    //           << "面板文件: " << panelFile << "\n"
    //           << "输出: " << outputFile << "\n"
    //           << "查询文件: " << (queryFile.empty() ? "未指定" : queryFile) << "\n"
    //           << "查询模式: " << queryMode << "\n"
    //           << "B: " << B << "\n"
    //           << "F: " << F << "\n"
    //           << "查询长度: " << queryLength << "\n";

    if (B==64)
    {
        // 执行程序逻辑
        wmFSPBWT<unsigned long long> CRY;
        CRY.F = F;
        CRY.B = B;
        CRY.T = pow(2, F);
        CRY.minSiteL = B * 2 - 1;

        int a = CRY.readMacsPanel(panelFile);
        std::cout << "读取面板: " << a << "\n";
        if (a != 0) return a;



        int c;
        if (queryMode == "in") {
            int b = CRY.makePanel();
            std::cout << "生成模糊面板: " << b << "\n";
            if (b != 0) return b;
            c = CRY.inPanelLongMatchQuery(queryLength, outputFile);
            std::cout << "面板内查询完成: " << c << "\n";
        } else if(queryMode == "out") {
            int w = CRY.readMacsQuery(queryFile);
            std::cout << "<读取查询完成>: " << w << "\n";
            int b = CRY.makePanel();
            std::cout << "生成模糊面板: " << b << "\n";
            if (b != 0) return b;
            c = CRY.outPanelLongMatchQuery(queryLength, outputFile);
            std::cout << "面板外查询完成: " << c << "\n";

        }
        string inf = panelFile +""+"_b"+ to_string(B)+"_f"+ to_string(F)+"_l"+ to_string(queryLength);
        CRY.outputInformationToFile(inf,queryMode);
    }
    else if (B==128)
    {
        // 执行程序逻辑
        wmFSPBWT<__uint128_t> CRY;
        CRY.F = F;
        CRY.B = B;
        CRY.T = pow(2, F);
        CRY.minSiteL = B * 2 - 1;

        int a = CRY.readMacsPanel(panelFile);
        std::cout << "读取面板: " << a << "\n";
        if (a != 0) return a;


        int c;
        if (queryMode == "in") {
            int b = CRY.makePanel();
            std::cout << "生成模糊面板: " << b << "\n";
            if (b != 0) return b;
            c = CRY.inPanelLongMatchQuery(queryLength, outputFile);
            std::cout << "面板内查询完成: " << c << "\n";
        } else if(queryMode == "out")
        {
            int w = CRY.readMacsQuery(queryFile);
            std::cout << "<读取查询完成>: " << w << "\n";
            int b = CRY.makePanel();
            std::cout << "生成模糊面板: " << b << "\n";
            if (b != 0) return b;
            c = CRY.outPanelLongMatchQuery(queryLength, outputFile);
            std::cout << "面板外查询完成: " << c << "\n";
        }
        string inf = panelFile +""+"_b"+ to_string(B)+"_f"+ to_string(F)+"_l"+ to_string(queryLength);
        CRY.outputInformationToFile(inf,queryMode);
    }
    return 0;
}