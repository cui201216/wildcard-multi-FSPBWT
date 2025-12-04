
#include <chrono>
#include<iostream>
#include <algorithm>
#include <array>
#include <cstring>
#include <fstream>
#include <numeric>
#include <sstream>
#include <vector>
#include<cmath>
#include <ctime>
#include <string>
#include"tools.h"
#include <unistd.h>
#include <unordered_map>
#include <utility>
#define FF 4

using namespace std;
namespace std {
    template<>
    struct hash<std::pair<int, int>> {
        size_t operator()(const std::pair<int, int>& p) const noexcept {
            // 结合两个 int 的哈希值
            size_t h1 = std::hash<int>{}(p.first);
            size_t h2 = std::hash<int>{}(p.second);
            // 简单组合哈希值（避免冲突）
            return h1 ^ (h2 + 0x9e3779b9 + (h1 << 6) + (h1 >> 2));
        }
    };
}
#define FF 4
using namespace std;
template<typename Syllable>
struct wmFSPBWT {
    int B;
    int F = 0;
    int T = 0;
    int M = 0;
    int N = 0;
    int n = 0;
    int minSiteL = 0;
    double readPanelTime = 0;
    double makeFuzzyPanelTime = 0;
    double inPanelQueryTime = 0;
    double inPanelFilteringTime = 0;
    double inPanelIdentificationTime = 0;
    double readQueryTime = 0;
    double outPanelQueryTime = 0;
    double outPanelFilteringTime = 0;
    double outPanelIdentificationTime = 0;
    u_long inPanelAlternativeNum = 0;
    u_long inPanelMatchNum = 0;
    u_long outPanelAlternativeNum = 0;
    u_long outPanelMatchNum = 0;
    u_long alternativeSyllableNum = 0;
    u_long matchLen = 0;
    vector<string> IDs;
    vector<int> physLocs;
    vector<vector<Syllable> > X; // MN bits
    vector<vector<uint32_t> > fuzzyX; // FMN/B bits
    vector<vector<int> > array; // 32MN/B bits
    vector<vector<int> > divergence; // 32MN/B bits

     vector<vector<bool>> panelMultiSyllable;
    vector<vector<bool>> queryMultiSyllable;
    int *u;

    int Q = 0;
    vector<vector<Syllable> > Z;
    vector<vector<uint32_t> > fuzzyZ;
    vector<string> qIDs;

    std::vector<std::vector<std::pair<unsigned int, uint8_t>>> panelMultiInfo; // M × n，存储 (start_index, length)
    std::vector<std::pair<uint8_t, uint8_t>> panelMultiValues; // 存储所有多字符位点 (position, value)

    std::vector<std::vector<std::pair<unsigned int, uint8_t>>> queryMultiInfo; // Q × n，存储 (start_index, length)
    std::vector<std::pair<uint8_t, uint8_t>> queryMultiValues; // 存储所有多字符位点 (position, value)

    u_long panelCount[10]={0};
    u_long queryCount[10]={0};


    vector<Syllable> filter;
    vector<vector<bool>> panelSyllableHavingMissing;
    vector<vector<bool>> querySyllableHavingMissing;
    std::unordered_map<std::pair<int,int>, Syllable> panelMissingData;
    std::unordered_map<std::pair<int,int>, Syllable> queryMissingData;
    int readMacsPanel(string txt_file);

    int readMacsQuery(string txt_file);

    int makePanel();

    int inPanelLongMatchQuery(int L, string inPanelOutput_file);

    int outPanelLongMatchQuery(int L, string outPanelOutput_file);

    bool inPanelSyllableMultiEqual(int index_a, int index_b, int k);

    void inPanelRefine(int L, int s_idx, int e_idx, int index_a, int index_b,
                       ofstream &out);

    void inPanelIdentification(int L, int s_idx, int e_idx, int index_a,
                               int index_b, ofstream &out);


    bool outPanelSyllableMultiEqual(int index_panel, int index_query, int k);

    void outPanelRefine(int L, int s_idx, int e_idx, int index_panel, int index_query,
                    ofstream &out);

    void outPanelIdentification(int L, int s_idx, int e_idx, int index_panel,
                                int index_query, ofstream &out);

    void outputInformationToFile(const std::string &fileName, string mode);

};


template<class Syllable>
int wmFSPBWT<Syllable>::readMacsPanel(string panel_file) {
    clock_t start, end;
    start = clock();
    std::ifstream in(panel_file);
    if (in.fail()) {
        std::cerr << "无法打开文件: " << panel_file << std::endl;
        return 1;
    }

    std::string line;

    // Step 1: 计算单倍型数 (M)
    M = 0;
    bool found_site = false;
    while (std::getline(in, line)) {
        if (line.rfind("SITE:", 0) == 0) {
            found_site = true;
            std::stringstream ss(line);
            std::string token;
            std::vector<std::string> tokens;
            while (std::getline(ss, token, '\t')) {
                tokens.push_back(token);
            }
            if (tokens.size() < 5) {
                std::cerr << "SITE行格式错误: 需要至少5个字段，实际为 " << tokens.size() << std::endl;
                return 2;
            }
            M = tokens[4].size();
            break;
        }
    }
    if (!found_site) {
        std::cerr << "未找到SITE行" << std::endl;
        return 2;
    }
    if (M < 1) {
        std::cerr << "无效的M: " << M << std::endl;
        return 3;
    }
    //std::cerr << "M = " << M << std::endl;

    // Step 2: 设置IDs
    IDs.resize(M);
    for (int i = 0; i < M; i++) {
        IDs[i] = std::to_string(i);
    }

    // Step 3: 计算位点数 (N)
    N = 0;
    in.clear();
    in.seekg(0);
    while (std::getline(in, line)) {
        if (line.rfind("SITE:", 0) == 0) {
            N++;
        }
    }
    if (N < 1) {
        std::cerr << "无效的N: " << N << std::endl;
        return 4;
    }
    n = (N + B - 1) / B;
    //std::cerr << "N = " << N << ", n = " << n << ", B = " << B << std::endl;

    // Step 4: 初始化数据结构
    try {
        X.resize(M, std::vector<Syllable>(n));
        panelMultiInfo.resize(M, std::vector<std::pair<unsigned int, uint8_t>>(n, std::make_pair(-1, 0)));
        //panelMultiValues.reserve(M*n*0.01);
        array.resize(n + 1, std::vector<int>(M));
        std::iota(array[0].begin(), array[0].end(), 0);
        divergence.resize(n + 1, std::vector<int>(M, 0));
        u = new int[(unsigned long)n * M * T];
        filter.resize(n);
        panelSyllableHavingMissing.resize(M,vector<bool>(n));
    } catch (const std::bad_alloc& e) {
        std::cerr << "内存分配失败: " << e.what() << std::endl;
        return -1;
    }

    // Step 5: 处理SITE行
    in.clear();
    in.seekg(0);
    std::vector<Syllable> X_(M, 0);
    std::vector<std::vector<std::pair<uint8_t, uint8_t>>> syllableMultis(M); // 临时存储多字符位点

    Syllable filterTemp=0;
    vector<Syllable> missingTemp;
    missingTemp.resize(M);
    Syllable one=1;

    int K = 0, k = 0;
    while (std::getline(in, line)) {
        if (line.rfind("SITE:", 0) != 0) {
            continue;
        }

        if (K >= N) {
            std::cerr << "SITE行数过多: K=" << K << ", 预期N=" << N << std::endl;
            return 10;
        }

        // 处理音节边界
        if (K % B == 0 && K != 0) {
            k = K / B - 1; // 上一个音节的索引
            if (k >= n) {
                std::cerr << "无效的k: " << k << ", n=" << n << std::endl;
                return 11;
            }
            for (int i = 0; i < M; i++) {
                X[i][k] = X_[i]; // 保存压缩音节
                // 存储多字符位点
                if (!syllableMultis[i].empty()) {
                    unsigned int start = panelMultiValues.size();
                    for (const auto& p : syllableMultis[i]) {
                        panelMultiValues.push_back(p);
                    }
                    panelMultiInfo[i][k] = std::make_pair(start, static_cast<uint8_t>(syllableMultis[i].size()));
                    syllableMultis[i].clear();
                } else {
                    panelMultiInfo[i][k] = std::make_pair(-1, 0);
                }
            }

            filter[k]=filterTemp;
            if (filterTemp!=0)
            {
                //having missing
                for (int i = 0; i < M; i++)
                {
                    if (panelSyllableHavingMissing[i][k]==true)
                    {
                        // add missing <i,k> , missingTemp[i]
                        panelMissingData[{i, k}] = missingTemp[i];
                    }
                }
            }
            X_.assign(M, 0); // 重置X_
            missingTemp.assign(M,0);
            filterTemp=0;
        }

        std::stringstream ss(line);
        std::string token;
        std::getline(ss, token, '\t'); // Skip "SITE:"
        std::getline(ss, token, '\t'); // Skip index
        std::getline(ss, token, '\t'); // Skip physLoc
        std::getline(ss, token, '\t'); // Skip other column
        std::getline(ss, token, '\t'); // Get haplotype data

        if (token.size() != M) {
            // ============== 新增调试代码开始 ==============
            std::cerr << "=========== DEBUG INFO ===========" << std::endl;
            std::cerr << "出错位置 (K): " << K << std::endl;
            std::cerr << "完整行内容: [" << line << "]" << std::endl;
            std::cerr << "提取的Token: [" << token << "]" << std::endl;
            std::cerr << "Token长度: " << token.size() << " (ASCI码: ";
            // 打印每个字符的ASCII码，防止有肉眼不可见的控制字符
            for(char c : token) {
                std::cerr << (int)c << " ";
            }
            std::cerr << ")" << std::endl;
            std::cerr << "==================================" << std::endl;
            // ============== 新增调试代码结束 ==============

            std::cerr << "单倍型数据长度不匹配: 预期 " << M << ", 实际 " << token.size() << ", K=" << K << std::endl;
            return 6;
        }

        int index = 0;
        for (char c : token) {
            if (index >= M) {
                std::cerr << "索引越界: index=" << index << ", M=" << M << ", K=" << K << std::endl;
                return 9;
            }
            if (K % B >= B) {
                std::cerr << "无效的K % B: " << K % B << ", K=" << K << std::endl;
                return 12;
            }
            if (c == '0') {
                X_[index] = X_[index] << 1;
            } else if (c == '1') {
                X_[index] = (X_[index] << 1) | 1;
            } else if (c > '1' && c <= '9') {
                X_[index] = (X_[index] << 1) | 1; // 多字符位点记为1
                syllableMultis[index].push_back(std::make_pair(static_cast<uint8_t>(K % B), static_cast<uint8_t>(c - '0')));
            } else if (c == '.') {
                X_[index] = (X_[index] << 1)|  1 ; // 缺失值位点记为1
                missingTemp[index]=missingTemp[index]|(one << ( B-1 - K%B) );
                panelSyllableHavingMissing[index][K/B]=true;
                filterTemp=filterTemp|(one<<( B-1 - K%B) );
                // // ---- 调试输出 ----
                // std::cerr << "[DEBUG] "
                //           << "i=" << index                // 样本编号
                //           << "  k=" << K/B              // 音节编号
                //           << "  K=" << K                  // 全局位点编号
                //           << "  K%B=" << K % B
                //           << "  bit-pos=" << (B - 1 - K % B)
                //           << std::endl;

            } else {
                std::cerr << "无效字符在 K=" << K << ", index=" << index << std::endl;
                return 7;
            }
            //panelCount[c - '0'] += 1;
            index++;
        }
        if (index != M) {
            std::cerr << "处理了 " << index << " 个单倍型，预期 " << M << ", K=" << K << std::endl;
            return 9;
        }

        // 处理最后一个音节
        if (K == N - 1) {
            k = K / B; // 最后一个音节的索引
            int pad2 = n * B - N; // 填充位点数
            if (pad2 < 0) {
                std::cerr << "无效的填充长度: pad2=" << pad2 << std::endl;
                return 14;
            }
            for (int i = 0; i < M; i++) {
                X_[i] <<= pad2; // 填充0
                X[i][k] = X_[i];
                if (!syllableMultis[i].empty()) {
                    unsigned int start = panelMultiValues.size();
                    for (const auto& p : syllableMultis[i]) {
                        panelMultiValues.push_back(p);
                    }
                    panelMultiInfo[i][k] = std::make_pair(start, static_cast<uint8_t>(syllableMultis[i].size()));
                    syllableMultis[i].clear();
                } else {
                    panelMultiInfo[i][k] = std::make_pair(-1, 0);
                }
                filter[k]=filterTemp;
                if (panelSyllableHavingMissing[i][k]==true) {
                    panelMissingData[{i, k}] = (missingTemp[i] << pad2);
                }
            }
        }
        K++;
    }

    if (K != N) {
        std::cerr << "处理了 " << K << " 个位点，预期 " << N << std::endl;
        return 10;
    }

    end = clock();
    readPanelTime = ((double)(end - start)) / CLOCKS_PER_SEC;
    std::cerr << "readPanelTime = " << readPanelTime << " 秒, 多字符位点数 = " << panelMultiValues.size() << std::endl;
    panelMultiValues.shrink_to_fit();
    return 0;
}

template<class Syllable>
int wmFSPBWT<Syllable>::readMacsQuery(string query_file) {
    clock_t start, end;
    start = clock();
    std::ifstream in(query_file);
    if (in.fail()) {
        std::cerr << "无法打开文件: " << query_file << std::endl;
        return 1;
    }

    std::string line;

    // Step 1: 计算单倍型数 (M)
    Q = 0;
    bool found_site = false;
    while (std::getline(in, line)) {
        if (line.rfind("SITE:", 0) == 0) {
            found_site = true;
            std::stringstream ss(line);
            std::string token;
            std::vector<std::string> tokens;
            while (std::getline(ss, token, '\t')) {
                tokens.push_back(token);
            }
            if (tokens.size() < 5) {
                std::cerr << "SITE行格式错误: 需要至少5个字段，实际为 " << tokens.size() << std::endl;
                return 2;
            }
            Q = tokens[4].size();
            break;
        }
    }
    if (!found_site) {
        std::cerr << "未找到SITE行" << std::endl;
        return 2;
    }
    if (Q < 1) {
        std::cerr << "无效的M: " << Q << std::endl;
        return 3;
    }
    std::cerr << "M = " << Q << std::endl;

    // Step 2: 设置IDs
    qIDs.resize(Q);
    for (int i = 0; i < Q; i++) {
        qIDs[i] = "query_" + std::to_string(i);
    }

    // Step 3: 计算位点数 (N)
    int siteNum = 0;
    in.clear();
    in.seekg(0);
    while (std::getline(in, line)) {
        if (line.rfind("SITE:", 0) == 0) {
            siteNum++;
        }
    }
    if (siteNum != N) {
        std::cerr << "panel跟query的位点数目不匹配 " << std::endl;
        return 4;
    }
    // Step 4: 初始化数据结构
    try {
        Z.resize(Q, std::vector<Syllable>(n));
        queryMultiInfo.resize(Q, std::vector<std::pair<unsigned int, uint8_t>>(n, std::make_pair(-1, 0)));
        //queryMultiValues.reserve(Q*n*0.01);
    } catch (const std::bad_alloc& e) {
        std::cerr << "内存分配失败: " << e.what() << std::endl;
        return -1;
    }

    // Step 5: 处理SITE行
    in.clear();
    in.seekg(0);
    std::vector<Syllable> Z_(M, 0);
    std::vector<std::vector<std::pair<uint8_t, uint8_t>>> syllableMultis(Q); // 临时存储多字符位点

    int K = 0, k = 0;
    while (std::getline(in, line)) {
        if (line.rfind("SITE:", 0) != 0) {
            continue;
        }

        if (K >= N) {
            std::cerr << "SITE行数过多: K=" << K << ", 预期N=" << N << std::endl;
            return 10;
        }

        // 处理音节边界
        if (K % B == 0 && K != 0) {
            k = K / B - 1; // 上一个音节的索引
            if (k >= n) {
                std::cerr << "无效的k: " << k << ", n=" << n << std::endl;
                return 11;
            }
            for (int i = 0; i < Q; i++) {
                Z[i][k] = Z_[i]; // 保存压缩音节
                // 存储多字符位点
                if (!syllableMultis[i].empty()) {
                    unsigned int start = queryMultiValues.size();
                    for (const auto& p : syllableMultis[i]) {
                        queryMultiValues.push_back(p);
                    }
                    queryMultiInfo[i][k] = std::make_pair(start, static_cast<uint8_t>(syllableMultis[i].size()));
                    syllableMultis[i].clear();
                } else {
                    queryMultiInfo[i][k] = std::make_pair(-1, 0);
                }
            }
            Z_.assign(M, 0); // 重置X_
        }

        std::stringstream ss(line);
        std::string token;
        std::getline(ss, token, '\t'); // Skip "SITE:"
        std::getline(ss, token, '\t'); // Skip index
        std::getline(ss, token, '\t'); // Skip physLoc
        std::getline(ss, token, '\t'); // Skip other column
        std::getline(ss, token, '\t'); // Get haplotype data

        if (token.size() != Q) {
            std::cerr << "单倍型数据长度不匹配: 预期 " << Q << ", 实际 " << token.size() << ", K=" << K << std::endl;
            return 6;
        }

        int index = 0;
        for (char c : token) {
            if (index >= Q) {
                std::cerr << "索引越界: index=" << index << ", M=" << M << ", K=" << K << std::endl;
                return 9;
            }
            if (c == '0') {
                Z_[index] = Z_[index] << 1;
            } else if (c == '1') {
                Z_[index] = (Z_[index] << 1) | 1;
            } else if (c > '1' && c <= '9') {
                Z_[index] = (Z_[index] << 1) | 1; // 多字符位点记为1
                syllableMultis[index].push_back(std::make_pair(static_cast<uint8_t>(K % B), static_cast<uint8_t>(c - '0')));
            } else if (c == '.') {
                std::cerr << "无效字符 '.' 在 K=" << K << ", index=" << index << std::endl;
                return 7;
            } else {
                std::cerr << "无效字符 '" << c << "' 在 K=" << K << ", index=" << index << std::endl;
                return 8;
            }
            queryCount[c - '0'] += 1;
            index++;
        }
        if (index != Q) {
            std::cerr << "处理了 " << index << " 个单倍型，预期 " << Q << ", K=" << K << std::endl;
            return 9;
        }

        // 处理最后一个音节
        if (K == N - 1) {
            k = K / B; // 最后一个音节的索引
            int pad2 = n * B - N; // 填充位点数

            for (int i = 0; i < Q; i++) {
                Z_[i] <<= pad2; // 填充0
                X[i][k] = Z_[i];
                if (!syllableMultis[i].empty()) {
                    unsigned int start = queryMultiValues.size();
                    for (const auto& p : syllableMultis[i]) {
                        queryMultiValues.push_back(p);
                    }
                    queryMultiInfo[i][k] = std::make_pair(start, static_cast<uint8_t>(syllableMultis[i].size()));
                    syllableMultis[i].clear();
                } else {
                    queryMultiInfo[i][k] = std::make_pair(-1, 0);
                }
            }
        }
        K++;
    }

    if (K != N) {
        std::cerr << "处理了 " << K << " 个位点，预期 " << N << std::endl;
        return 10;
    }
    queryMultiValues.shrink_to_fit();
    end = clock();
    readQueryTime = ((double)(end - start)) / CLOCKS_PER_SEC;
    std::cerr << "readPanelTime = " << readPanelTime << " 秒, 多字符位点数 = " << panelMultiValues.size() << std::endl;
    return 0;
}

template<class Syllable>
int wmFSPBWT<Syllable>::makePanel() // fuzzy way :overall
{
    //get fuzzyX
    clock_t start, end;
    start = clock();
    if (FF * n % 32 != 0) {
        fuzzyX.resize(M, vector<uint32_t>(FF * n / 32 + 1));
    } else {
        fuzzyX.resize(M, vector<uint32_t>(FF * n / 32));
    }

    // // 输出 fuzzyX 的大小
    // std::cout << "fuzzyX 的行数（M）：" << fuzzyX.size() << std::endl;
    // if (!fuzzyX.empty()) {
    //     std::cout << "fuzzyX 的列数（每个内部 vector 的大小）：" << fuzzyX[0].size() << std::endl;
    // }
    //
    // // 计算 fuzzyX 的总内存占用
    // size_t total_size = 0;
    // for (const auto& row : fuzzyX) {
    //     total_size += row.size() * sizeof(uint32_t);
    // }
    // std::cout << "fuzzyX 的总内存占用（字节）：" << total_size << std::endl;


    fuzzyX.shrink_to_fit();
    //
    // // 输出 fuzzyX 的大小
    // std::cout << "fuzzyX 的行数（M）：" << fuzzyX.size() << std::endl;
    // if (!fuzzyX.empty()) {
    //     std::cout << "fuzzyX 的列数（每个内部 vector 的大小）：" << fuzzyX[0].size() << std::endl;
    // }
    //
    // // 计算 fuzzyX 的总内存占用
    // total_size = 0;
    // for (const auto& row : fuzzyX) {
    //     total_size += row.size() * sizeof(uint32_t);
    // }
    // std::cout << "fuzzyX 的总内存占用（字节）：" << total_size << std::endl;

    int index = 0, count = 0;
    //int fuzzyIndex = 0;
    for (int k = 0; k < n; k++) {
        if (count == 32 / FF) {
            ++index;
            count = 0;
        }
        for (int i = 0; i < M; i++) {
            int num = 0;
            if (B == 64) {
                num = __builtin_popcountll(X[i][k]);
            } else if (B == 128) {
                num = countSetBits128(X[i][k]);
            }
            uint32_t temp = num % T;
            fuzzyX[i][index] = (fuzzyX[i][index] << FF) | temp;
        }
        count++;
    }
    int pad1 = 32 - count * FF;
    for (int i = 0; i < M; i++) {
        fuzzyX[i][fuzzyX[0].size() - 1] = fuzzyX[i][fuzzyX[0].size() - 1] << pad1;
    }

    //get array[][] divergence[][] u[][][]
    //u.resize(n, vector<vector<int>>(M, vector<int>(T)));
    int a_count[T] = {0};
    int d_count[T] = {0};
    int a[T][M];
    int d[T][M];
    //vector<int> a[T]; // make array[][]  his size for u[][]
    //vector<int> d[T]; //  make divergence
    //vector<int> p(T);	//  make divergence
    int p[T];
    uint32_t temp, fuzzy;
    int unit_index;
    int offset;
    int index1;
    int m;
    int plus;
    for (int k = 0; k < n; k++) //each position
    {
        for (int _ = 0; _ < T; _++) {
            p[_] = k + 1;
        }
        unit_index = k / (32 / FF);
        offset = 32 - FF - k % (32 / FF) * FF; //get fuzzy from fuzzyX

        for (int i = 0; i < M; i++) //each hapolotype
        {
            //update u[][][]
            for (int _ = 0; _ < T; _++) {
                u[k * (M * T) + i * T + _]=a_count[_];
                //u[k][i][_] = a_count[_];
                if (divergence[k][i] > p[_]) {
                    p[_] = divergence[k][i];
                }
            }
            index1 = array[k][i];
            //uint32_t temp = fuzzyX[index][unit_index];

            temp = fuzzyX[index1][unit_index];
            if (F == 1) {
                fuzzy = (temp >> offset) & 0b1;
            } else if (F == 2) {
                fuzzy = (temp >> offset) & 0b11;
            } else if (F == 3) {
                fuzzy = (temp >> offset) & 0b111;
            } else if (F == 4) {
                fuzzy = (temp >> offset) & 0b1111;
            }

            a[fuzzy][a_count[fuzzy]] = index1;
            d[fuzzy][d_count[fuzzy]] = p[fuzzy];
            a_count[fuzzy]++;
            d_count[fuzzy]++;
            p[fuzzy] = 0;
        } //end j in  M
        //update array divergence u[][][]
        m = 0;
        for (int _ = 0; _ < T; _++) {
            for (int w = 0; w < a_count[_]; w++) {
                array[k + 1][m] = a[_][w];
                divergence[k + 1][m] = d[_][w];
                ++m;
            }
        }
        if (m != M) {
            return 4;
        }
        //put plus into u for w()
        for (int _ = 1; _ < T; _++) {
            plus = 0;
            for (int j = 0; j < _; j++) {
                plus += a_count[j];
            }
            for (int j = 0; j < M; j++) {
                u[k * (M * T) + j * T + _] += plus;
            }
        }
        for (int _ = 0; _ < T; _++) {
            a_count[_] = 0;
            d_count[_] = 0;
        }
    }
    end = clock();
    makeFuzzyPanelTime = ((double) (end - start)) / CLOCKS_PER_SEC;
    return 0;
}

template<class Syllable>
int wmFSPBWT<Syllable>::inPanelLongMatchQuery(int L, string inPanelOutput_file) {
    clock_t start, end;
    start = clock();
    if (L < minSiteL) {
        return 1;
    }
    ofstream out(inPanelOutput_file);
    if (out.fail())
        return 2;

    int l = (L - (B - 1)) / B;
    int k;
    for (k = 0; k < n - 1; k++) {
        int unit_index = k / (32 / FF);
        int offset = 32 - FF - k % (32 / FF) * FF;

        bool m[T];
        for (int _ = 0; _ < T; _++) {
            m[_] = false;
        }
        int top = 0;
        bool report = false;
        for (int i = 0; i < M; i++) {
            if (divergence[k][i] > k - l) {
                for (int w = 0; w < T - 1; w++) {
                    for (int v = w + 1; v < T; v++) {
                        if (m[w] == true && m[v] == true) {
                            report = true;
                            break;
                        }
                    }
                }
                if (report == true) {
                    for (int i_a = top; i_a < i - 1; i_a++) {
                        int maxDivergence = 0;
                        for (int i_b = i_a + 1; i_b < i; i_b++) {
                            if (divergence[k][i_b] > maxDivergence) {
                                maxDivergence = divergence[k][i_b];
                            }
                            uint32_t temp1, temp2, fuzzy1, fuzzy2;
                            int index_a = array[k][i_a], index_b =
                                    array[k][i_b];
                            temp1 = fuzzyX[index_a][unit_index];
                            temp2 = fuzzyX[index_b][unit_index];
                            if (F == 1) {
                                fuzzy1 = (temp1 >> offset) & 0b1;
                                fuzzy2 = (temp2 >> offset) & 0b1;
                            } else if (F == 2) {
                                fuzzy1 = (temp1 >> offset) & 0b11;
                                fuzzy2 = (temp2 >> offset) & 0b11;
                            } else if (F == 3) {
                                fuzzy1 = (temp1 >> offset) & 0b111;
                                fuzzy2 = (temp2 >> offset) & 0b111;
                            } else if (F == 4) {
                                fuzzy1 = (temp1 >> offset) & 0b1111;
                                fuzzy2 = (temp2 >> offset) & 0b1111;
                            }

                            if (fuzzy1 != fuzzy2) {
                                inPanelIdentification(L, maxDivergence - 1,
                                                      k, index_a, index_b, out);
                                ++inPanelAlternativeNum;
                            }
                        }
                    }
                    report = false;
                } //end if(report==true)

                top = i;
                for (int _ = 0; _ < T; _++) {
                    m[_] = false;
                }
                //report = false;
            } // end if(divergence[k][i]>k-l)
            //change m[]

            /*int fuzzy;
             fuzzy = catchFuzzy(array[k][i], k);*/

            uint32_t temp, fuzzy;
            temp = fuzzyX[array[k][i]][unit_index];
            if (F == 1) {
                fuzzy = (temp >> offset) & 0b1;
            } else if (F == 2) {
                fuzzy = (temp >> offset) & 0b11;
            } else if (F == 3) {
                fuzzy = (temp >> offset) & 0b111;
            } else if (F == 4) {
                fuzzy = (temp >> offset) & 0b1111;
            }

            m[fuzzy] = true;
        } //for i from 0 to M-1
        //cheak bottom block
        for (int w = 0; w < T - 1; w++) {
            for (int v = w + 1; v < T; v++) {
                if (m[w] == true && m[v] == true) {
                    report = true;
                }
            }
        }
        if (report == true) {
            for (int i_a = top; i_a < M - 1; i_a++) {
                int maxDivergence = 0;
                for (int i_b = i_a + 1; i_b < M; i_b++) {
                    if (divergence[k][i_b] > maxDivergence) {
                        maxDivergence = divergence[k][i_b];
                    }
                    uint32_t temp1, temp2, fuzzy1, fuzzy2;
                    int index_a = array[k][i_a];
                    int index_b = array[k][i_b];
                    temp1 = fuzzyX[index_a][unit_index];
                    temp2 = fuzzyX[index_b][unit_index];
                    if (F == 1) {
                        fuzzy1 = (temp1 >> offset) & 0b1;
                        fuzzy2 = (temp2 >> offset) & 0b1;
                    } else if (F == 2) {
                        fuzzy1 = (temp1 >> offset) & 0b11;
                        fuzzy2 = (temp2 >> offset) & 0b11;
                    } else if (F == 3) {
                        fuzzy1 = (temp1 >> offset) & 0b111;
                        fuzzy2 = (temp2 >> offset) & 0b111;
                    } else if (F == 4) {
                        fuzzy1 = (temp1 >> offset) & 0b1111;
                        fuzzy2 = (temp2 >> offset) & 0b1111;
                    }

                    if (fuzzy1 != fuzzy2) {
                        inPanelIdentification(L, maxDivergence - 1, k,
                                              index_a, index_b, out);
                        ++inPanelAlternativeNum;
                    }
                }
            }
        }
    }
    int top = 0;
    for (int i = 0; i < M; i++) {
        if (divergence[k][i] > k - l + 1) {
            for (int i_a = top; i_a < i - 1; i_a++) {
                int maxDivergence = 0;
                for (int i_b = i_a + 1; i_b < i; i_b++) {
                    if (divergence[k][i_b] > maxDivergence) {
                        maxDivergence = divergence[k][i_b];
                    }
                    int unit_index = k / (32 / FF);
                    int offset = 32 - FF - k % (32 / FF) * FF;
                    //int fuzzy1, fuzzy2;
                    int index_a = array[k][i_a];
                    int index_b = array[k][i_b];
                    uint32_t temp1, temp2, fuzzy1, fuzzy2;
                    //fuzzy1 = catchFuzzy(index_a, k);
                    //fuzzy2 = catchFuzzy(index_b, k);
                    temp1 = fuzzyX[index_a][unit_index];
                    temp2 = fuzzyX[index_b][unit_index];

                    if (F == 1) {
                        fuzzy1 = (temp1 >> offset) & 0b1;
                        fuzzy2 = (temp2 >> offset) & 0b1;
                    } else if (F == 2) {
                        fuzzy1 = (temp1 >> offset) & 0b11;
                        fuzzy2 = (temp2 >> offset) & 0b11;
                    } else if (F == 3) {
                        fuzzy1 = (temp1 >> offset) & 0b111;
                        fuzzy2 = (temp2 >> offset) & 0b111;
                    } else if (F == 4) {
                        fuzzy1 = (temp1 >> offset) & 0b1111;
                        fuzzy2 = (temp2 >> offset) & 0b1111;
                    }

                    if (fuzzy1 == fuzzy2) {
                        inPanelIdentification(L, maxDivergence - 1, n,
                                              index_a, index_b, out);
                        ++inPanelAlternativeNum;
                    } else if (fuzzy1 != fuzzy2) {
                        if (k - maxDivergence >= l) {
                            inPanelIdentification(L, maxDivergence - 1,
                                                  n - 1, index_a, index_b, out);
                            ++inPanelAlternativeNum;
                        }
                    }
                }
            }
            top = i;
        }
    }
    //bottom block
    for (int i_a = top; i_a < M - 1; i_a++) {
        int maxDivergence = 0;
        for (int i_b = i_a + 1; i_b < M; i_b++) {
            int index_a = array[k][i_a];
            int index_b = array[k][i_b];
            if (divergence[k][i_b] > maxDivergence) {
                maxDivergence = divergence[k][i_b];
            }
            inPanelIdentification(L, maxDivergence - 1, n, index_a, index_b,
                                  out);
            ++inPanelAlternativeNum;
        }
    }

    end = clock();

    inPanelQueryTime = ((double) (end - start)) / CLOCKS_PER_SEC;

    inPanelFilteringTime = inPanelQueryTime - inPanelIdentificationTime;
    return 0;
    out.close();
    cout << "matches has been put into " << inPanelOutput_file << endl;
}

template<class Syllable>
int wmFSPBWT<Syllable>::outPanelLongMatchQuery(int L, string outPanelOutput_file) {
    clock_t start, end;
    start = clock();

    if (FF * n % 32 != 0) {
        fuzzyZ.resize(Q, vector<uint32_t>(FF * n / 32 + 1));
    } else {
        fuzzyZ.resize(Q, vector<uint32_t>(FF * n / 32));
    }
        int index = 0, count = 0;
        for (int k = 0; k < n; k++) {
            if (count == 32 / FF) {
                ++index;
                count = 0;
            }
            for (int i = 0; i < Q; i++) {
                int num = 0;
                if (B == 64) {
                    num = __builtin_popcountll(Z[i][k]);
                } else if (B == 128) {
                    num = countSetBits128(Z[i][k]);
                }
                uint32_t temp = num % T;
                fuzzyZ[i][index] = (fuzzyZ[i][index] << FF) | temp;
            }
            count++;
        }
        int pad1 = 32 - count * FF;
        for (int i = 0; i < Q; i++) {
            fuzzyZ[i][fuzzyZ[0].size() - 1] =
                    fuzzyZ[i][fuzzyZ[0].size() - 1] << pad1;
        }

    if (L < minSiteL) {
        return 1;
    }
    int l = (L - (B - 1)) / B;
    ofstream out(outPanelOutput_file);
    if (out.fail())
        return 2;

    vector<int> dZ(M); //match start
    dZ.shrink_to_fit();
    vector<int> t(n + 1); //fake location
    t.shrink_to_fit();
    vector<int> Zdivergence(n + 2); //divergence of Z	因为要从n到0计算Zdivergence和belowZdivergence时要用到dZ[n+1]=n
    Zdivergence.shrink_to_fit();
    vector<int> belowZdivergence(n + 2); //divergence of Z.below
    belowZdivergence.shrink_to_fit();
    for (int q = 0; q < Q; q++) {
        fill(dZ.begin(), dZ.end(), 0);
        fill(t.begin(), t.end(), 0);
        fill(Zdivergence.begin(), Zdivergence.end(), 0);
        fill(belowZdivergence.begin(), belowZdivergence.end(), 0);

        string qID = qIDs[q >> 1] + "-" + to_string(q & 1);
        t[0] = 0;

        //fake location
        for (int k = 0; k < n; k++) {
            uint32_t temp, fuzzy;
            int unit_index = k / (32 / FF);
            int offset = 32 - FF - k % (32 / FF) * FF;
            temp = fuzzyZ[q][unit_index];

            if (F == 1) {
                fuzzy = (temp >> offset) & 0b1;
            } else if (F == 2) {
                fuzzy = (temp >> offset) & 0b11;
            } else if (F == 3) {
                fuzzy = (temp >> offset) & 0b111;
            } else if (F == 4) {
                fuzzy = (temp >> offset) & 0b1111;
            }

            if (t[k] != M) {
                t[k + 1] =
                    u[k * (M * T) + t[k] * T + fuzzy];
//                    u[k][t[k]][fuzzy];
            } else {
                // t[k] == M
                if (fuzzy < T - 1) {
                    t[k + 1] =u[k * (M * T) + 0 * T + fuzzy+1];
                        //u[k][0][fuzzy + 1];
                } else if (fuzzy == T - 1) {
                    t[k + 1] = M;
                } else {
                    return 2;
                }
            }
        }

        Zdivergence[n + 1] = belowZdivergence[n + 1] = n;
        for (int k = n; k >= 0; --k) {
            Zdivergence[k] = std::min(Zdivergence[k + 1], k);
            belowZdivergence[k] = std::min(belowZdivergence[k + 1], k);
            if (t[k] != 0) {
                uint32_t panelTemp, panelFuzzy, queryTemp, queryFuzzy;
                int unit_index = (Zdivergence[k] - 1) / (32 / FF);
                int reminder = (Zdivergence[k] - 1) % (32 / FF);
                int offset = 32 - FF - reminder * FF;

                int index = array[k][t[k] - 1];
                //int index = array[Zdivergence[k] - 1][t[k] - 1];	//hapolotype on query
                panelTemp = fuzzyX[index][unit_index];
                queryTemp = fuzzyZ[q][unit_index];
                if (F == 1) {
                    panelFuzzy = (panelTemp >> offset) & 0b1;
                    queryFuzzy = (queryTemp >> offset) & 0b1;
                } else if (F == 2) {
                    panelFuzzy = (panelTemp >> offset) & 0b11;
                    queryFuzzy = (queryTemp >> offset) & 0b11;
                } else if (F == 3) {
                    panelFuzzy = (panelTemp >> offset) & 0b111;
                    queryFuzzy = (queryTemp >> offset) & 0b111;
                } else if (F == 4) {
                    panelFuzzy = (panelTemp >> offset) & 0b1111;
                    queryFuzzy = (queryTemp >> offset) & 0b1111;
                }
                //向前更新Zdivergence
                while (Zdivergence[k] > 0 && panelFuzzy == queryFuzzy) {
                    --Zdivergence[k];
                    if (reminder == 0) {
                        if (unit_index == 0) {
                            break;
                        } else {
                            reminder = 32 / FF - 1; //last bug is here ^^
                            --unit_index;
                        }
                    } else {
                        --reminder;
                    }
                    offset = 32 - FF - reminder * FF;
                    panelTemp = fuzzyX[index][unit_index];
                    queryTemp = fuzzyZ[q][unit_index];
                    if (F == 1) {
                        panelFuzzy = (panelTemp >> offset) & 0b1;
                        queryFuzzy = (queryTemp >> offset) & 0b1;
                    } else if (F == 2) {
                        panelFuzzy = (panelTemp >> offset) & 0b11;
                        queryFuzzy = (queryTemp >> offset) & 0b11;
                    } else if (F == 3) {
                        panelFuzzy = (panelTemp >> offset) & 0b111;
                        queryFuzzy = (queryTemp >> offset) & 0b111;
                    } else if (F == 4) {
                        panelFuzzy = (panelTemp >> offset) & 0b1111;
                        queryFuzzy = (queryTemp >> offset) & 0b1111;
                    }
                }
            } else {
                //t[k]==0
                Zdivergence[k] = k;
            }
            if (t[k] < M) {
                uint32_t panelTemp, panelFuzzy, queryTemp, queryFuzzy;
                int unit_index = (belowZdivergence[k] - 1) / (32 / FF);
                int reminder = (belowZdivergence[k] - 1) % (32 / FF);
                int offset = 32 - FF - reminder * FF;
                int index = array[k][t[k]]; //hapolotype below query
                panelTemp = fuzzyX[index][unit_index];
                queryTemp = fuzzyZ[q][unit_index];
                if (F == 1) {
                    panelFuzzy = (panelTemp >> offset) & 0b1;
                    queryFuzzy = (queryTemp >> offset) & 0b1;
                } else if (F == 2) {
                    panelFuzzy = (panelTemp >> offset) & 0b11;
                    queryFuzzy = (queryTemp >> offset) & 0b11;
                } else if (F == 3) {
                    panelFuzzy = (panelTemp >> offset) & 0b111;
                    queryFuzzy = (queryTemp >> offset) & 0b111;
                } else if (F == 4) {
                    panelFuzzy = (panelTemp >> offset) & 0b1111;
                    queryFuzzy = (queryTemp >> offset) & 0b1111;
                }
                //向前更新belowZdivergence

                while (belowZdivergence[k] > 0 && panelFuzzy == queryFuzzy) {
                    belowZdivergence[k]--;
                    if (reminder == 0) {
                        if (unit_index == 0) {
                            break;
                        }
                        reminder = 32 / FF - 1;
                        --unit_index;
                    } else {
                        --reminder;
                    }

                    offset = 32 - FF - reminder * FF;
                    panelTemp = fuzzyX[index][unit_index];
                    queryTemp = fuzzyZ[q][unit_index];
                    if (F == 1) {
                        panelFuzzy = (panelTemp >> offset) & 0b1;
                        queryFuzzy = (queryTemp >> offset) & 0b1;
                    } else if (F == 2) {
                        panelFuzzy = (panelTemp >> offset) & 0b11;
                        queryFuzzy = (queryTemp >> offset) & 0b11;
                    } else if (F == 3) {
                        panelFuzzy = (panelTemp >> offset) & 0b111;
                        queryFuzzy = (queryTemp >> offset) & 0b111;
                    } else if (F == 4) {
                        panelFuzzy = (panelTemp >> offset) & 0b1111;
                        queryFuzzy = (queryTemp >> offset) & 0b1111;
                    }
                }
            } else {
                belowZdivergence[k] = k;
            }
        }

        int f, g;
        f = g = t[0];
        vector<int> ftemp, gtemp;
        ftemp.resize(T);
        gtemp.resize(T);

        for (int k = 0; k < n; k++) {
            int unit_index = k / (32 / FF);
            int offset = 32 - FF - k % (32 / FF) * FF;
            uint32_t temp = fuzzyZ[q][unit_index];
            uint32_t fuzzyQ;
            if (F == 1) {
                fuzzyQ = (temp >> offset) & 0b1;
            } else if (F == 2) {
                fuzzyQ = (temp >> offset) & 0b11;
            } else if (F == 3) {
                fuzzyQ = (temp >> offset) & 0b111;
            } else if (F == 4) {
                fuzzyQ = (temp >> offset) & 0b1111;
            }
            if (g == M) {
                if (f == M) {
                    //update ftemp
                    for (int i = 0; i < T; i++) {
                        if (fuzzyQ != i) {
                            if (i != T - 1) {
                                ftemp[i] = u[k * (M * T) + 0 * T + i+1];
                                   // u[k][0][i + 1];
                            } else {
                                ftemp[i] = M;
                            }
                        }
                    }

                    if (fuzzyQ != T - 1) {
                        f =u[k * (M * T) + 0 * T + fuzzyQ+1];
                            //u[k][0][fuzzyQ + 1];
                    } else {
                        f = M;
                    }
                } else //f!=M
                {
                    for (int i = 0; i < T; i++) {
                        if (fuzzyQ != i) {
                            ftemp[i] =u[k * (M * T) + f * T + i];
                                //u[k][f][i];
                        }
                    }
                    f = u[k * (M * T) + f * T + fuzzyQ];
                        //u[k][f][fuzzyQ];
                }
                //update gtemp and g
                for (int i = 0; i < T; i++) {
                    if (fuzzyQ != i) {
                        if (i < T - 1) {
                            gtemp[i] =u[k * (M * T) + 0 * T + i + 1];
                                //u[k][0][i + 1];
                        } else {
                            gtemp[i] = M;
                        }
                    }
                }
                if (fuzzyQ < T - 1) {
                    g =u[k * (M * T) + 0 * T + fuzzyQ+1];
                        //u[k][0][fuzzyQ + 1];
                } else {
                    g = M;
                }
            } else //g!=M
            {
                for (int i = 0; i < T; i++) {
                    if (i != fuzzyQ) {
                        ftemp[i] = u[k * (M * T) + f * T + i];
                            //u[k][f][i];
                        gtemp[i] = u[k * (M * T) + g * T + i];
                            //u[k][g][i];
                    }
                }
                f =u[k * (M * T) + f * T + fuzzyQ];
                    //u[k][f][fuzzyQ];
                g =u[k * (M * T) + g * T + fuzzyQ];
                    //u[k][g][fuzzyQ];
            }

            //output matches
            for (int i = 0; i < T; i++) {
                if (i != fuzzyQ) {
                    while (ftemp[i] != gtemp[i]) {
                        //output Match
                        //int start = 0, end = 0;
                        int index = array[k + 1][ftemp[i]];
                        outPanelIdentification(L, dZ[index] - 1, k, index,
                                               q, out);
                        ++outPanelAlternativeNum;
                        ++ftemp[i];
                    }
                }
            }

            if (f == g) {
                if (k + 1 - Zdivergence[k + 1] == l) {
                    --f;
                    dZ[array[k + 1][f]] = k + 1 - l;
                    //store divergence
                }

                //if (k + 1 - belowZdivergence[k + 1] == l)
                if (k + 1 - belowZdivergence[k + 1] == l) {
                    //store divergence
                    dZ[array[k + 1][g]] = k + 1 - l;
                    ++g;
                }
            }
            if (f != g) {
                while (divergence[k + 1][f] <= k + 1 - l) {
                    --f;
                    dZ[array[k + 1][f]] = k + 1 - l;
                }
                while (g < M && divergence[k + 1][g] <= k + 1 - l) {
                    dZ[array[k + 1][g]] = k + 1 - l;

                    ++g;
                }
            }
        }

        //mathces no ending at
        while (f != g) {
            //output Match
            //	int start = 0, end = 0;
            int index = array[n][f];

            outPanelIdentification(L, dZ[index] - 1, n, index, q, out);
            ++outPanelAlternativeNum;

            ++f;
        }

        std::fill(t.begin(), t.end(), 0);
        std::fill(Zdivergence.begin(), Zdivergence.end(), 0);
        std::fill(belowZdivergence.begin(), belowZdivergence.end(), 0);
    }

    end = clock();

    outPanelQueryTime = ((double) (end - start)) / CLOCKS_PER_SEC;

    outPanelFilteringTime = outPanelQueryTime - outPanelIdentificationTime;
    cout << "matches has been put into " << outPanelOutput_file << endl;
    return 0;
}

template<class Syllable>
bool wmFSPBWT<Syllable>::inPanelSyllableMultiEqual(int index_a, int index_b, int k) {
    if (X[index_a][k] != X[index_b][k]) {
        return false;
    }
    auto info_a = panelMultiInfo[index_a][k];
    auto info_b = panelMultiInfo[index_b][k];
    if (info_a.second != info_b.second) {
        return false;
    }
    if (info_a.second == 0) {
        return true;
    }
    // 比较多字符位点列表
    for (unsigned int i = 0; i < info_a.second; i++) {
        if (panelMultiValues[info_a.first + i] != panelMultiValues[info_b.first + i]) {
            return false;
        }
    }
    return true;
}

template<class Syllable>
bool wmFSPBWT<Syllable>::outPanelSyllableMultiEqual(int index_panel, int index_query, int k) {
    if (X[index_panel][k] != Z[index_query][k]) {
        return false;
    }
    auto info_panel = panelMultiInfo[index_panel][k];
    auto info_query = queryMultiInfo[index_query][k];
    if (info_panel.second != info_query.second) {
        return false;
    }
    if (info_panel.second == 0) {
        return true;
    }
    // 比较多字符位点列表
    for (unsigned int i = 0; i < info_panel.second; i++) {
        if (panelMultiValues[info_panel.first + i] != queryMultiValues[info_query.first + i]) {
            return false;
        }
    }
    return true;
}


template<class Syllable>
void wmFSPBWT<Syllable>::inPanelRefine(int L, int s_idx, int e_idx, int index_a, int index_b,
                                          std::ofstream &out) {
    bool startDone=false;
    bool endDone=false;
    int start = 0, end = 0;
    if (s_idx == -1) {
        start = 0;
        startDone = true;
    } else {
        if (panelMultiInfo[index_a][s_idx].second == 0 && panelMultiInfo[index_b][s_idx].second == 0) {
            unsigned long tz;
            if (B == 64) {
                tz = __builtin_ctzll(X[index_a][s_idx] ^ X[index_b][s_idx]);
            } else if (B == 128) {
                tz = ctz128_uint128(X[index_a][s_idx] ^ X[index_b][s_idx]);
            }
            start = (s_idx + 1) * B - tz;
            startDone = true;
        } else if (panelMultiInfo[index_a][s_idx].second > 0 && panelMultiInfo[index_b][s_idx].second > 0) {
            auto info_a = panelMultiInfo[index_a][s_idx];
            auto info_b = panelMultiInfo[index_b][s_idx];
            if (info_a.second == info_b.second) {
                bool equal = true;
                for (unsigned int i = 0; i < info_a.second; i++) {
                    if (panelMultiValues[info_a.first + i] != panelMultiValues[info_b.first + i]) {
                        equal = false;
                        break;
                    }
                }
                if (equal) {
                    unsigned long tz;
                    if (B == 64) {
                        tz = __builtin_ctzll(X[index_a][s_idx] ^ X[index_b][s_idx]);
                    } else if (B == 128) {
                        tz = ctz128_uint128(X[index_a][s_idx] ^ X[index_b][s_idx]);
                    }
                    start = (s_idx + 1) * B - tz;
                    startDone = true;
                }
            }
        }
        if (startDone==false)
        {
            auto info_a = panelMultiInfo[index_a][s_idx];
            auto info_b = panelMultiInfo[index_b][s_idx];
            // 构建音节的位点序列
            std::vector<char> seq_a(B, '0'), seq_b(B, '0');
            // 填充二进制位点
            for (int j = 0; j < B; ++j) {
                seq_a[j] = ((X[index_a][s_idx] >> (B - 1 - j)) & 1) ? '1' : '0';
                seq_b[j] = ((X[index_b][s_idx] >> (B - 1 - j)) & 1) ? '1' : '0';
            }
            // 填充多字符位点
            if (info_a.second > 0) {
                for (unsigned int i = 0; i < info_a.second; ++i) {
                    int pos = panelMultiValues[info_a.first + i].first;
                    char val = '0' + panelMultiValues[info_a.first + i].second;
                    seq_a[pos] = val;
                }
            }
            if (info_b.second > 0) {
                for (unsigned int i = 0; i < info_b.second; ++i) {
                    int pos = panelMultiValues[info_b.first + i].first;
                    char val = '0' + panelMultiValues[info_b.first + i].second;
                    seq_b[pos] = val;
                }
            }
            // 从音节末尾向前比较
            int tz = 0;
            for (int j = B - 1; j >= 0; --j) {
                if (seq_a[j] != seq_b[j]) {
                    break;
                }
                ++tz;
            }
            start = (s_idx + 1) * B - tz;
            startDone = true;
        }

    }

    if (e_idx == n) {
        end = N;
        endDone = true;
    } else {
        if (panelMultiInfo[index_a][e_idx].second == 0 && panelMultiInfo[index_b][e_idx].second == 0) {
            unsigned long tz = 0;
            if (B == 64) {
                tz = __builtin_clzll(X[index_a][e_idx] ^ X[index_b][e_idx]);
            } else if (B == 128) {
                tz = clz128_uint128(X[index_a][e_idx] ^ X[index_b][e_idx]);
            }
            end = e_idx * B + tz;
            endDone=true;
        } else if (panelMultiInfo[index_a][e_idx].second > 0 && panelMultiInfo[index_b][e_idx].second > 0) {
            auto info_a = panelMultiInfo[index_a][e_idx];
            auto info_b = panelMultiInfo[index_b][e_idx];
            if (info_a.second == info_b.second) {
                bool equal = true;
                for (unsigned int i = 0; i < info_a.second; i++) {
                    if (panelMultiValues[info_a.first + i] != panelMultiValues[info_b.first + i]) {
                        equal = false;
                        break;
                    }
                }
                if (equal) {
                    unsigned long tz = 0;
                    if (B == 64) {
                        tz = __builtin_clzll(X[index_a][e_idx] ^ X[index_b][e_idx]);
                    } else if (B == 128) {
                        tz = clz128_uint128(X[index_a][e_idx] ^ X[index_b][e_idx]);
                    }
                    end = e_idx * B + tz;
                    endDone=true;
                }
            }
        }
        if (endDone==false)
        {
            auto info_a = panelMultiInfo[index_a][e_idx];
            auto info_b = panelMultiInfo[index_b][e_idx];
            std::vector<char> seq_a(B, '0'), seq_b(B, '0');
            for (int j = 0; j < B; ++j) {
                seq_a[j] = ((X[index_a][e_idx] >> (B - 1 - j)) & 1) ? '1' : '0';
                seq_b[j] = ((X[index_b][e_idx] >> (B - 1 - j)) & 1) ? '1' : '0';
            }
            if (info_a.second > 0) {
                for (unsigned int i = 0; i < info_a.second; ++i) {
                    int pos = panelMultiValues[info_a.first + i].first;
                    char val = '0' + panelMultiValues[info_a.first + i].second;
                    seq_a[pos] = val;
                }
            }
            if (info_b.second > 0) {
                for (unsigned int i = 0; i < info_b.second; ++i) {
                    int pos = panelMultiValues[info_b.first + i].first;
                    char val = '0' + panelMultiValues[info_b.first + i].second;
                    seq_b[pos] = val;
                }
            }
            // 从音节开头向后比较
            int prefix = 0;
            for (int j = 0; j < B; ++j) {
                if (seq_a[j] != seq_b[j]) {
                    break;
                }
                ++prefix;
            }
            end = e_idx * B + prefix;
            endDone=true;
        }
    }

    if (end - start >= L) {
        out << IDs[index_a] << '\t' << IDs[index_b] << '\t' << start << '\t'
            << end-1 << '\n';
        ++inPanelMatchNum;
        matchLen += (end - start);
    }
}

template<class Syllable>
void wmFSPBWT<Syllable>::outPanelRefine(int L, int s_idx, int e_idx, int index_panel, int index_query,
                                          std::ofstream &out) {
    bool startDone=false;
    bool endDone=false;
    int start = 0, end = 0;
    if (s_idx == -1) {
        start = 0;
        startDone = true;
    } else {
        if (panelMultiInfo[index_panel][s_idx].second == 0 && queryMultiInfo[index_query][s_idx].second == 0) {
            unsigned long tz;
            if (B == 64) {
                tz = __builtin_ctzll(X[index_panel][s_idx] ^ Z[index_query][s_idx]);
            } else if (B == 128) {
                tz = ctz128_uint128(X[index_panel][s_idx] ^ Z[index_query][s_idx]);
            }
            start = (s_idx + 1) * B - tz;
            startDone = true;
        } else if (panelMultiInfo[index_panel][s_idx].second > 0 && queryMultiInfo[index_query][s_idx].second > 0) {
            auto info_a = panelMultiInfo[index_panel][s_idx];
            auto info_b = queryMultiInfo[index_query][s_idx];
            if (info_a.second == info_b.second) {
                bool equal = true;
                for (unsigned int i = 0; i < info_a.second; i++) {
                    if (panelMultiValues[info_a.first + i] != queryMultiValues[info_b.first + i]) {
                        equal = false;
                        break;
                    }
                }
                if (equal) {
                    unsigned long tz;
                    if (B == 64) {
                        tz = __builtin_ctzll(X[index_panel][s_idx] ^ Z[index_query][s_idx]);
                    } else if (B == 128) {
                        tz = ctz128_uint128(X[index_panel][s_idx] ^ Z[index_query][s_idx]);
                    }
                    start = (s_idx + 1) * B - tz;
                    startDone = true;
                }
            }
        }
        if (startDone==false)
        {
            auto info_a = panelMultiInfo[index_panel][s_idx];
            auto info_b = queryMultiInfo[index_query][s_idx];
            // 构建音节的位点序列
            std::vector<char> seq_a(B, '0'), seq_b(B, '0');
            // 填充二进制位点
            for (int j = 0; j < B; ++j) {
                seq_a[j] = ((X[index_panel][s_idx] >> (B - 1 - j)) & 1) ? '1' : '0';
                seq_b[j] = ((Z[index_query][s_idx] >> (B - 1 - j)) & 1) ? '1' : '0';
            }
            // 填充多字符位点
            if (info_a.second > 0) {
                for (unsigned int i = 0; i < info_a.second; ++i) {
                    int pos = panelMultiValues[info_a.first + i].first;
                    char val = '0' + panelMultiValues[info_a.first + i].second;
                    seq_a[pos] = val;
                }
            }
            if (info_b.second > 0) {
                for (unsigned int i = 0; i < info_b.second; ++i) {
                    int pos = queryMultiValues[info_b.first + i].first;
                    char val = '0' + queryMultiValues[info_b.first + i].second;
                    seq_b[pos] = val;
                }
            }
            // 从音节末尾向前比较
            int tz = 0;
            for (int j = B - 1; j >= 0; --j) {
                if (seq_a[j] != seq_b[j]) {
                    break;
                }
                ++tz;
            }
            start = (s_idx + 1) * B - tz;
            startDone = true;
        }

    }

    if (e_idx == n) {
        end = N;
        endDone = true;
    } else {
        if (panelMultiInfo[index_panel][e_idx].second == 0 && queryMultiInfo[index_query][e_idx].second == 0) {
            unsigned long tz = 0;
            if (B == 64) {
                tz = __builtin_clzll(X[index_panel][e_idx] ^ Z[index_query][e_idx]);
            } else if (B == 128) {
                tz = clz128_uint128(X[index_panel][e_idx] ^ Z[index_query][e_idx]);
            }
            end = e_idx * B + tz;
            endDone=true;
        } else if (panelMultiInfo[index_panel][e_idx].second > 0 && queryMultiInfo[index_query][e_idx].second > 0) {
            auto info_a = panelMultiInfo[index_panel][e_idx];
            auto info_b = queryMultiInfo[index_query][e_idx];
            if (info_a.second == info_b.second) {
                bool equal = true;
                for (unsigned int i = 0; i < info_a.second; i++) {
                    if (panelMultiValues[info_a.first + i] != queryMultiValues[info_b.first + i]) {
                        equal = false;
                        break;
                    }
                }
                if (equal) {
                    unsigned long tz = 0;
                    if (B == 64) {
                        tz = __builtin_clzll(X[index_panel][e_idx] ^ Z[index_query][e_idx]);
                    } else if (B == 128) {
                        tz = clz128_uint128(X[index_panel][e_idx] ^ Z[index_query][e_idx]);
                    }
                    end = e_idx * B + tz;
                    endDone=true;
                }
            }
        }
        if (endDone==false)
        {
            auto info_a = panelMultiInfo[index_panel][e_idx];
            auto info_b = queryMultiInfo[index_query][e_idx];
            std::vector<char> seq_a(B, '0'), seq_b(B, '0');
            for (int j = 0; j < B; ++j) {
                seq_a[j] = ((X[index_panel][e_idx] >> (B - 1 - j)) & 1) ? '1' : '0';
                seq_b[j] = ((Z[index_query][e_idx] >> (B - 1 - j)) & 1) ? '1' : '0';
            }
            if (info_a.second > 0) {
                for (unsigned int i = 0; i < info_a.second; ++i) {
                    int pos = panelMultiValues[info_a.first + i].first;
                    char val = '0' + panelMultiValues[info_a.first + i].second;
                    seq_a[pos] = val;
                }
            }
            if (info_b.second > 0) {
                for (unsigned int i = 0; i < info_b.second; ++i) {
                    int pos = queryMultiValues[info_b.first + i].first;
                    char val = '0' + queryMultiValues[info_b.first + i].second;
                    seq_b[pos] = val;
                }
            }
            // 从音节开头向后比较
            int prefix = 0;
            for (int j = 0; j < B; ++j) {
                if (seq_a[j] != seq_b[j]) {
                    break;
                }
                ++prefix;
            }
            end = e_idx * B + prefix;
            endDone=true;
        }
    }

    if (end - start >= L) {
        out << IDs[index_panel] << '\t' << qIDs[index_query] << '\t' << start << '\t'
            << end-1 << '\n';
        ++outPanelMatchNum;
        matchLen += (end - start);
    }
}


template<class Syllable>
void wmFSPBWT<Syllable>::inPanelIdentification(int L, int s_idx, int e_idx, int index_a,
                                                int index_b, std::ofstream &out) {
    alternativeSyllableNum += (e_idx - s_idx + 1);
    clock_t start, end;
    start = clock();
    int l = (L - 2 * (B - 1)) / B;
    int head = s_idx + 1, tail;
    while (head < e_idx) {
        tail = head;
        while (tail < e_idx && inPanelSyllableMultiEqual(index_a, index_b, tail)) {
            tail++;
        }
        if (tail - head >= l) {
            inPanelRefine(L, head - 1, tail, index_a, index_b, out);
        }
        head = tail + 1;
        while (head < e_idx && !inPanelSyllableMultiEqual(index_a, index_b, head)) {
            head++;
        }
    }
    end = clock();
    inPanelIdentificationTime += ((double) (end - start)) / CLOCKS_PER_SEC;
}

template<class Syllable>
void wmFSPBWT<Syllable>::outPanelIdentification(int L, int s_idx, int e_idx, int index_panel,
                                                int index_query, std::ofstream &out) {
    alternativeSyllableNum += (e_idx - s_idx + 1);
    clock_t start, end;
    start = clock();
    int l = (L - 2 * (B - 1)) / B;
    int head = s_idx + 1, tail;
    while (head < e_idx) {
        tail = head;
        while (tail < e_idx && outPanelSyllableMultiEqual(index_panel, index_query, tail)) {
            tail++;
        }
        if (tail - head >= l) {
            outPanelRefine(L, head - 1, tail, index_panel, index_query, out);
        }
        head = tail + 1;
        while (head < e_idx && !outPanelSyllableMultiEqual(index_panel, index_query, head)) {
            head++;
        }
    }
    end = clock();
    outPanelIdentificationTime += ((double) (end - start)) / CLOCKS_PER_SEC;
}

template<class Syllable>
void wmFSPBWT<Syllable>::outputInformationToFile(const std::string &fileName, string mode) {
    std::ofstream outputFile(fileName);

    if (outputFile.is_open()) {
        if (mode == "in") {
            outputFile << "B: " << B << std::endl;
            outputFile << "F: " << F << std::endl;
            outputFile << "T: " << T << std::endl;
            outputFile << "M: " << M << std::endl;
            outputFile << "N: " << N << std::endl;
            outputFile << "n: " << n << std::endl;
            outputFile << "minSiteL: " << minSiteL << std::endl;
            outputFile << endl;

            outputFile << "readPanelTime: " << readPanelTime << std::endl;
            outputFile << "makeFuzzyPanelTime: " << makeFuzzyPanelTime
                    << std::endl;
            outputFile << "inPanelQueryTime: " << inPanelQueryTime
                    << std::endl;
            outputFile << "inPanelFilteringTime: " << inPanelFilteringTime
                    << std::endl;
            outputFile << "inPanelIdentificationTime: "
                    << inPanelIdentificationTime << std::endl;
            outputFile << "inPanelAlternativeNum: " << inPanelAlternativeNum
                    << std::endl;
            outputFile << "inPanelMatchNum: " << inPanelMatchNum
                    << std::endl;
            outputFile << "AlternativeLen: " << alternativeSyllableNum*B
        << std::endl;
            outputFile << "MatchLen: " << matchLen
                    << std::endl;
            // 新增：输出面板位点计数信息
            outputFile << "Panel Site Counts:" << std::endl;
            for (int i = 0; i < 10; i++) {
                outputFile << "Site " << i << ": " << panelCount[i] << std::endl;
            }
            std::cout << "Information has been written to " << fileName
                    << std::endl;
        } else if (mode == "out") {
            outputFile << "B: " << B << std::endl;
            outputFile << "F: " << F << std::endl;
            outputFile << "T: " << T << std::endl;
            outputFile << "M: " << M << std::endl;
            outputFile << "N: " << N << std::endl;
            outputFile << "n: " << n << std::endl;
            outputFile << "Q: " << Q << std::endl;
            outputFile << "minSiteL: " << minSiteL << std::endl;
            outputFile << endl;

            outputFile << "readPanelTime: " << readPanelTime << std::endl;
            outputFile << "makeFuzzyPanelTime: " << makeFuzzyPanelTime
                    << std::endl;

            outputFile << "readQueryTime: " << readQueryTime << std::endl;
            outputFile << "outPanelQueryTime: " << outPanelQueryTime
                    << std::endl;
            outputFile << "outPanelFilteringTime: " << outPanelFilteringTime
                    << std::endl;
            outputFile << "outPanelIdentificationTime: "
                    << outPanelIdentificationTime << std::endl;
            outputFile << endl;

            outputFile << "outPanelAlternativeNum: "
                    << outPanelAlternativeNum << std::endl;
            outputFile << "outPanelMatchNum: " << outPanelMatchNum
                    << std::endl;
            outputFile << "AlternativeLen: " << alternativeSyllableNum*B
<< std::endl;
            outputFile << "MatchLen: " << matchLen
                    << std::endl;
            // 新增：输出面板和查询位点计数信息
            outputFile << "Panel Site Counts:" << std::endl;
            for (int i = 0; i < 10; i++) {
                outputFile << "Site " << i << ": " << panelCount[i] << std::endl;
            }
            outputFile << std::endl;

            outputFile << "Query Site Counts:" << std::endl;
            for (int i = 0; i < 10; i++) {
                outputFile << "Site " << i << ": " << queryCount[i] << std::endl;
            }
            outputFile << std::endl;
            outputFile.close();
            std::cout << "Information has been written to " << fileName
                    << std::endl;
        }
    } else {
        std::cerr << "Unable to open file " << fileName << " for writing."
                << std::endl;
    }
}

