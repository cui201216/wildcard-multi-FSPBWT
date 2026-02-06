#include <chrono>
#include<iostream>

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

namespace std
{
    template <>
    struct hash<std::pair<int, int>>
    {
        size_t operator()(const std::pair<int, int>& p) const noexcept
        {
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

template <typename Syllable>
struct wmFSPBWT
{
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
    vector<vector<Syllable>> X; // MN bits
    vector<vector<uint32_t>> fuzzyX; // FMN/B bits
    vector<vector<int>> array; // 32MN/B bits
    vector<vector<int>> divergence; // 32MN/B bits

    vector<vector<bool>> panelMultiSyllable;
    vector<vector<bool>> queryMultiSyllable;
    int* u;

    int Q = 0;
    vector<vector<Syllable>> Z;
    vector<vector<uint32_t>> fuzzyZ;
    vector<string> qIDs;

    std::vector<std::vector<std::pair<unsigned int, uint8_t>>> panelMultiInfo; // M × n，存储 (start_index, length)
    std::vector<std::pair<uint8_t, uint8_t>> panelMultiValues; // 存储所有多字符位点 (position, value)

    std::vector<std::vector<std::pair<unsigned int, uint8_t>>> queryMultiInfo; // Q × n，存储 (start_index, length)
    std::vector<std::pair<uint8_t, uint8_t>> queryMultiValues; // 存储所有多字符位点 (position, value)

    u_long panelCount[11] = {0};
    u_long queryCount[11] = {0};


    vector<Syllable> filter;
    vector<vector<bool>> panelSyllableHavingMissing;
    vector<vector<bool>> querySyllableHavingMissing;
    std::unordered_map<std::pair<int, int>, Syllable> panelMissingData;
    std::unordered_map<std::pair<int, int>, Syllable> queryMissingData;

    int readVcfPanel(string panel_file);

    int readVcfQuery(string query_file);

    int readMacsPanel(string txt_file);

    int readMacsQuery(string txt_file);

    int makePanel();

    int inPanelLongMatchQuery(int L, string inPanelOutput_file);

    int outPanelLongMatchQuery(int L, string outPanelOutput_file);

    bool inPanelSyllableEqual(int index_a, int index_b, int k);

    void inPanelRefine(int L, int s_idx, int e_idx, int index_a, int index_b,
                       ofstream& out);

    void inPanelIdentification(int L, int s_idx, int e_idx, int index_a,
                               int index_b, ofstream& out);


    bool outPanelSyllableEqual(int index_panel, int index_query, int k);

    void outPanelRefine(int L, int s_idx, int e_idx, int index_panel, int index_query,
                        ofstream& out);

    void outPanelIdentification(int L, int s_idx, int e_idx, int index_panel,
                                int index_query, ofstream& out);

    void outputInformationToFile(const std::string& fileName, string mode);

    void save(const string& index_file);
    void load(const string& index_file);
};

template <class Syllable>
int wmFSPBWT<Syllable>::readVcfPanel(string panel_file)
{
    clock_t start, end;
    start = clock();

    // 使用 C 风格 IO 提速
    FILE* fp = fopen(panel_file.c_str(), "r");
    if (!fp) return 1;

    char* buffer = new char[1024 * 1024]; // 1MB 行缓存，足够大
    size_t buf_size = 1024 * 1024;

    // 1. 跳过 Meta headers (##)
    long data_start_pos = 0;
    while (getline(&buffer, &buf_size, fp) != -1) {
        if (buffer[0] != '#' || buffer[1] != '#') {
            break; // 读到了 #CHROM 行
        }
    }

    // 2. 解析 Header (#CHROM) 计算样本数
    // 逻辑：计算第9个tab后面的列数
    int tab_count = 0;
    char* ptr = buffer;
    while (*ptr) {
        if (*ptr == '\t') tab_count++;
        ptr++;
    }
    int sample_count = tab_count - 8; // 前9列是Header，后面是样本
    if (sample_count < 1) return 2;
    M = sample_count * 2; // VCF 默认双单倍型

    // 3. 计算行数 (N)
    data_start_pos = ftell(fp); // 记录数据开始位置
    N = 0;
    // 快速过一遍文件数行数
    while (getline(&buffer, &buf_size, fp) != -1) {
        if (buffer[0] != '#') N++; // 简单判断
    }
    if (N < 1) return 2;

    n = (N + B - 1) / B;

    // 4. 初始化结构
    try {
        IDs.resize(M);
        for(int i=0; i<M; ++i) IDs[i] = to_string(i); // 简化ID
        physLocs.resize(N);
        X.resize(M, vector<Syllable>(n));

        panelMultiInfo.resize(M, vector<pair<unsigned int, uint8_t>>(n, make_pair(-1, 0)));
        array.resize(n + 1, vector<int>(M));
        divergence.resize(n + 1, vector<int>(M, 0));
        u = new int[(long)n * M * T];
        filter.resize(n, 0);
        panelSyllableHavingMissing.resize(M, vector<bool>(n, false));

        iota(array[0].begin(), array[0].end(), 0);
    } catch (...) { return -1; }

    // 重置文件指针
    fseek(fp, data_start_pos, SEEK_SET);

    vector<Syllable> X_(M, 0);
    vector<vector<pair<uint8_t, uint8_t>>> syllableMultis(M);
    vector<Syllable> missingTemp(M, 0);
    Syllable filterTemp = 0;
    Syllable one = 1;

    // 5. 极速读取循环
    for (int K = 0, k = 0; K < N; K++)
    {
        if (getline(&buffer, &buf_size, fp) == -1) break;

        k = K / B;

        // --- Block 保存逻辑 (保持不变) ---
        if (K % B == 0 && K != 0) {
            int prev_k = k - 1;
            for (int i = 0; i < M; i++) {
                X[i][prev_k] = X_[i];
                if (!syllableMultis[i].empty()) {
                    unsigned int start_idx = panelMultiValues.size();
                    for (const auto& p : syllableMultis[i]) panelMultiValues.push_back(p);
                    panelMultiInfo[i][prev_k] = make_pair(start_idx, (uint8_t)syllableMultis[i].size());
                    syllableMultis[i].clear();
                } else {
                    panelMultiInfo[i][prev_k] = make_pair(-1, 0);
                }
                if (panelSyllableHavingMissing[i][prev_k]) {
                    panelMissingData[{i, prev_k}] = missingTemp[i];
                }
            }
            filter[prev_k] |= filterTemp;
            fill(X_.begin(), X_.end(), 0);
            fill(missingTemp.begin(), missingTemp.end(), 0);
            filterTemp = 0;
        }

        // --- 核心解析逻辑：跳过前9列，然后只认字符 ---
        char* p = buffer;
        int tabs = 0;

        // 快速跳过前9列，顺便抓取 POS (第2列，索引1)
        while (tabs < 9 && *p) {
            if (tabs == 1) { // POS列
                physLocs[K] = atoi(p); // 读数字
                // 移动指针直到下一个 tab
                while (*p && *p != '\t') p++;
                if (*p == '\t') { tabs++; p++; continue; }
            }
            if (*p == '\t') tabs++;
            p++;
        }

        // 扫描基因型
        int index = 0;
        while (*p && *p != '\n' && index < M) {
            char c = *p;

            if (c >= '0' && c <= '9') {
                // 处理数字位点
                int val = c - '0';
                if (val == 0) {
                    X_[index] = X_[index] << 1;
                    panelCount[0]++;
                } else if (val == 1) {
                    X_[index] = (X_[index] << 1) | 1;
                    panelCount[1]++;
                } else {
                    // val > 1 (2-9)
                    X_[index] = (X_[index] << 1) | 1; // 模糊化为1
                    panelCount[val]++;
                    syllableMultis[index].push_back(make_pair((uint8_t)(K % B), (uint8_t)val));
                }
                index++;
            }
            else if (c == '.') {
                // 处理缺失值
                X_[index] = (X_[index] << 1) | 1; // 模糊化为1
                panelCount[10]++;
                missingTemp[index] |= (one << (B - 1 - K % B));
                panelSyllableHavingMissing[index][k] = true;
                filterTemp |= (one << (B - 1 - K % B));
                index++;
            }
            // 其他字符 (| / \t :) 全部忽略
            p++;
        }

        // --- 最后一个 Block 填充 (保持不变) ---
        if (K == N - 1) {
            int pad2 = n * B - N;
            for (int i = 0; i < M; i++) {
                if (pad2 > 0) X_[i] <<= pad2;
                X[i][k] = X_[i];
                if (!syllableMultis[i].empty()) {
                    unsigned int start_idx = panelMultiValues.size();
                    for (const auto& kv : syllableMultis[i]) panelMultiValues.push_back(kv);
                    panelMultiInfo[i][k] = make_pair(start_idx, (uint8_t)syllableMultis[i].size());
                } else {
                    panelMultiInfo[i][k] = make_pair(-1, 0);
                }
                if (panelSyllableHavingMissing[i][k]) {
                    panelMissingData[{i, k}] = missingTemp[i];
                }
            }
            filter[k] |= filterTemp;
        }
    }

    free(buffer);
    fclose(fp);
    panelMultiValues.shrink_to_fit();
    end = clock();
    readPanelTime = ((double)(end - start)) / CLOCKS_PER_SEC;
    return 0;
}

template <class Syllable>
int wmFSPBWT<Syllable>::readVcfQuery(string query_file)
{
    clock_t start, end;
    start = clock();

    FILE* fp = fopen(query_file.c_str(), "r");
    if (!fp) return 1;

    char* buffer = new char[1024 * 1024];
    size_t buf_size = 1024 * 1024;

    // 1. Skip Headers
    while (getline(&buffer, &buf_size, fp) != -1) {
        if (buffer[0] != '#' || buffer[1] != '#') break;
    }

    // 2. Count Samples
    int tab_count = 0;
    char* ptr = buffer;
    while (*ptr) {
        if (*ptr == '\t') tab_count++;
        ptr++;
    }
    int sample_count = tab_count - 8;
    if (sample_count < 1) return 2;
    Q = sample_count * 2;

    // 3. Init
    try {
        qIDs.resize(Q);
        for(int i=0; i<Q; ++i) qIDs[i] = "Q" + to_string(i);
        Z.resize(Q, vector<Syllable>(n));
        queryMultiInfo.resize(Q, vector<pair<unsigned int, uint8_t>>(n, make_pair(-1, 0)));
        querySyllableHavingMissing.resize(Q, vector<bool>(n, false));
    } catch (...) { return -1; }

    long data_pos = ftell(fp); // 数据开始处

    vector<Syllable> Z_(Q, 0);
    vector<vector<pair<uint8_t, uint8_t>>> syllableMultis(Q);
    vector<Syllable> missingTemp(Q, 0);
    Syllable one = 1;

    // 4. 读取数据
    for (int K = 0, k = 0; K < N; K++)
    {
        if (getline(&buffer, &buf_size, fp) == -1) break;
        k = K / B;

        // Block处理
        if (K % B == 0 && K != 0) {
            int prev_k = k - 1;
            for (int i = 0; i < Q; i++) {
                Z[i][prev_k] = Z_[i];
                if (!syllableMultis[i].empty()) {
                    unsigned int start_idx = queryMultiValues.size();
                    for (const auto& p : syllableMultis[i]) queryMultiValues.push_back(p);
                    queryMultiInfo[i][prev_k] = make_pair(start_idx, (uint8_t)syllableMultis[i].size());
                    syllableMultis[i].clear();
                } else {
                    queryMultiInfo[i][prev_k] = make_pair(-1, 0);
                }
                if (missingTemp[i] != 0) {
                    queryMissingData[{i, prev_k}] = missingTemp[i];
                    filter[prev_k] |= missingTemp[i];
                }
            }
            fill(Z_.begin(), Z_.end(), 0);
            fill(missingTemp.begin(), missingTemp.end(), 0);
        }

        // 解析逻辑：跳过前9列，扫描字符
        char* p = buffer;
        int tabs = 0;
        // 快速跳9个Tab
        while (tabs < 9 && *p) {
            if (*p == '\t') tabs++;
            p++;
        }

        int index = 0;
        while (*p && *p != '\n' && index < Q) {
            char c = *p;
            if (c >= '0' && c <= '9') {
                int val = c - '0';
                if (val == 0) {
                    Z_[index] = Z_[index] << 1;
                    queryCount[0]++;
                } else if (val == 1) {
                    Z_[index] = (Z_[index] << 1) | 1;
                    queryCount[1]++;
                } else { // 2-9
                    Z_[index] = (Z_[index] << 1) | 1;
                    queryCount[val]++;
                    syllableMultis[index].push_back(make_pair((uint8_t)(K % B), (uint8_t)val));
                }
                index++;
            } else if (c == '.') {
                Z_[index] = (Z_[index] << 1) | 1;
                queryCount[10]++;
                missingTemp[index] |= (one << (B - 1 - K % B));
                querySyllableHavingMissing[index][k] = true;
                index++;
            }
            p++;
        }

        // 最后一个Block
        if (K == N - 1) {
            int pad2 = n * B - N;
            for (int i = 0; i < Q; i++) {
                if (pad2 > 0) Z_[i] <<= pad2;
                Z[i][k] = Z_[i];
                if (!syllableMultis[i].empty()) {
                    unsigned int start_idx = queryMultiValues.size();
                    for (const auto& kv : syllableMultis[i]) queryMultiValues.push_back(kv);
                    queryMultiInfo[i][k] = make_pair(start_idx, (uint8_t)syllableMultis[i].size());
                } else {
                    queryMultiInfo[i][k] = make_pair(-1, 0);
                }
                if (querySyllableHavingMissing[i][k]) {
                    queryMissingData[{i, k}] = missingTemp[i];
                    filter[k] |= missingTemp[i];
                }
            }
        }
    }

    free(buffer);
    fclose(fp);
    queryMultiValues.shrink_to_fit();
    end = clock();
    readQueryTime = ((double)(end - start)) / CLOCKS_PER_SEC;
    return 0;
}
template <class Syllable>
int wmFSPBWT<Syllable>::readMacsPanel(string panel_file)
{
    clock_t start, end;
    start = clock();
    std::ifstream in(panel_file);
    if (in.fail())
    {
        std::cerr << "无法打开文件: " << panel_file << std::endl;
        return 1;
    }

    std::string line;

    // Step 1: 计算单倍型数 (M)
    M = 0;
    bool found_site = false;
    while (std::getline(in, line))
    {
        if (line.rfind("SITE:", 0) == 0)
        {
            found_site = true;
            std::stringstream ss(line);
            std::string token;
            std::vector<std::string> tokens;
            while (std::getline(ss, token, '\t'))
            {
                tokens.push_back(token);
            }
            if (tokens.size() < 5)
            {
                std::cerr << "SITE行格式错误: 需要至少5个字段，实际为 " << tokens.size() << std::endl;
                return 2;
            }
            M = tokens[4].size();
            break;
        }
    }
    if (!found_site)
    {
        std::cerr << "未找到SITE行" << std::endl;
        return 2;
    }
    if (M < 1)
    {
        std::cerr << "无效的M: " << M << std::endl;
        return 3;
    }
    //std::cerr << "M = " << M << std::endl;

    // Step 2: 设置IDs
    IDs.resize(M);
    for (int i = 0; i < M; i++)
    {
        IDs[i] = std::to_string(i);
    }

    // Step 3: 计算位点数 (N)
    N = 0;
    in.clear();
    in.seekg(0);
    while (std::getline(in, line))
    {
        if (line.rfind("SITE:", 0) == 0)
        {
            N++;
        }
    }
    if (N < 1)
    {
        std::cerr << "无效的N: " << N << std::endl;
        return 4;
    }
    n = (N + B - 1) / B;
    //std::cerr << "N = " << N << ", n = " << n << ", B = " << B << std::endl;

    // Step 4: 初始化数据结构
    try
    {
        X.resize(M, std::vector<Syllable>(n));
        panelMultiInfo.resize(M, std::vector<std::pair<unsigned int, uint8_t>>(n, std::make_pair(-1, 0)));
        //panelMultiValues.reserve(M*n*0.01);
        array.resize(n + 1, std::vector<int>(M));
        std::iota(array[0].begin(), array[0].end(), 0);
        divergence.resize(n + 1, std::vector<int>(M, 0));
        u = new int[(unsigned long)n * M * T];
        filter.resize(n);
        panelSyllableHavingMissing.resize(M, vector<bool>(n));
    }
    catch (const std::bad_alloc& e)
    {
        std::cerr << "内存分配失败: " << e.what() << std::endl;
        return -1;
    }

    // Step 5: 处理SITE行
    in.clear();
    in.seekg(0);
    std::vector<Syllable> X_(M, 0);
    std::vector<std::vector<std::pair<uint8_t, uint8_t>>> syllableMultis(M); // 临时存储多字符位点

    Syllable filterTemp = 0;
    vector<Syllable> missingTemp;
    missingTemp.resize(M);
    Syllable one = 1;

    int K = 0, k = 0;
    while (std::getline(in, line))
    {
        if (line.rfind("SITE:", 0) != 0)
        {
            continue;
        }

        if (K >= N)
        {
            std::cerr << "SITE行数过多: K=" << K << ", 预期N=" << N << std::endl;
            return 10;
        }

        // 处理音节边界
        if (K % B == 0 && K != 0)
        {

            k = K / B - 1; // 上一个音节的索引
            if (k >= n)
            {
                std::cerr << "无效的k: " << k << ", n=" << n << std::endl;
                return 11;
            }
            for (int i = 0; i < M; i++)
            {
                X[i][k] = X_[i]; // 保存压缩音节
                // 存储多字符位点
                if (!syllableMultis[i].empty())
                {
                    unsigned int start = panelMultiValues.size();
                    for (const auto& p : syllableMultis[i])
                    {
                        panelMultiValues.push_back(p);
                    }
                    panelMultiInfo[i][k] = std::make_pair(start, static_cast<uint8_t>(syllableMultis[i].size()));
                    syllableMultis[i].clear();
                }
                else
                {
                    panelMultiInfo[i][k] = std::make_pair(-1, 0);
                }
            }

            filter[k] = filterTemp;
            if (filterTemp != 0)
            {
                //having missing
                for (int i = 0; i < M; i++)
                {
                    if (panelSyllableHavingMissing[i][k] == true)
                    {
                        // add missing <i,k> , missingTemp[i]
                        panelMissingData[{i, k}] = missingTemp[i];
                    }
                }
            }
            X_.assign(M, 0); // 重置X_
            missingTemp.assign(M, 0);
            filterTemp = 0;

        }

        std::stringstream ss(line);
        std::string token;
        std::getline(ss, token, '\t'); // Skip "SITE:"
        std::getline(ss, token, '\t'); // Skip index
        std::getline(ss, token, '\t'); // Skip physLoc
        std::getline(ss, token, '\t'); // Skip other column
        std::getline(ss, token, '\t'); // Get haplotype data

        if (token.size() != M)
        {
            std::cerr << "单倍型数据长度不匹配: 预期 " << M << ", 实际 " << token.size() << ", K=" << K << std::endl;
            return 6;
        }

        int index = 0;
        for (char c : token)
        {
            if (index >= M)
            {
                std::cerr << "索引越界: index=" << index << ", M=" << M << ", K=" << K << std::endl;
                return 9;
            }
            if (K % B >= B)
            {
                std::cerr << "无效的K % B: " << K % B << ", K=" << K << std::endl;
                return 12;
            }
            if (c == '0')
            {
                X_[index] = X_[index] << 1;
                ++panelCount[0];
            }
            else if (c == '1')
            {
                X_[index] = (X_[index] << 1) | 1;
                ++panelCount[1];
            }
            else if (c > '1' && c <= '9')
            {
                X_[index] = (X_[index] << 1) | 1; // 多字符位点记为1
                ++panelCount[c-'0'];
                syllableMultis[index].push_back(
                    std::make_pair(static_cast<uint8_t>(K % B), static_cast<uint8_t>(c - '0')));
            }
            else if (c == '.')
            {
                X_[index] = (X_[index] << 1) | 1; // 缺失值位点记为1
                ++panelCount[10];
                missingTemp[index] = missingTemp[index] | (one << (B - 1 - K % B));
                panelSyllableHavingMissing[index][K / B] = true;
                filterTemp = filterTemp | (one << (B - 1 - K % B));
            }
            else
            {
                std::cerr << "无效字符在 K=" << K << ", index=" << index << std::endl;
                return 7;
            }
            //panelCount[c - '0'] += 1;
            index++;
        }
        if (index != M)
        {
            std::cerr << "处理了 " << index << " 个单倍型，预期 " << M << ", K=" << K << std::endl;
            return 9;
        }

        // 处理最后一个音节
        if (K == N - 1)
        {
            k = K / B; // 最后一个音节的索引
            int pad2 = n * B - N; // 填充位点数
            if (pad2 < 0)
            {
                std::cerr << "无效的填充长度: pad2=" << pad2 << std::endl;
                return 14;
            }
            for (int i = 0; i < M; i++)
            {
                X_[i] <<= pad2; // 填充0
                X[i][k] = X_[i];
                if (!syllableMultis[i].empty())
                {
                    unsigned int start = panelMultiValues.size();
                    for (const auto& p : syllableMultis[i])
                    {
                        panelMultiValues.push_back(p);
                    }
                    panelMultiInfo[i][k] = std::make_pair(start, static_cast<uint8_t>(syllableMultis[i].size()));
                    syllableMultis[i].clear();
                }
                else
                {
                    panelMultiInfo[i][k] = std::make_pair(-1, 0);
                }
                filter[k] = filterTemp;
                if (panelSyllableHavingMissing[i][k] == true)
                {
                    panelMissingData[{i, k}] = missingTemp[i];
                }
            }
        }
        K++;
    }

    if (K != N)
    {
        std::cerr << "处理了 " << K << " 个位点，预期 " << N << std::endl;
        return 10;
    }

    end = clock();
    readPanelTime = ((double)(end - start)) / CLOCKS_PER_SEC;
    panelMultiValues.shrink_to_fit();
    return 0;
}

template <class Syllable>
int wmFSPBWT<Syllable>::readMacsQuery(string query_file)
{
    clock_t start, end;
    start = clock();
    std::ifstream in(query_file);
    if (in.fail())
    {
        std::cerr << "无法打开文件: " << query_file << std::endl;
        return 1;
    }

    std::string line;

    // Step 1: 计算单倍型数 (Q)
    Q = 0;
    bool found_site = false;
    while (std::getline(in, line))
    {
        if (line.rfind("SITE:", 0) == 0)
        {
            found_site = true;
            std::stringstream ss(line);
            std::string token;
            std::vector<std::string> tokens;
            while (std::getline(ss, token, '\t'))
            {
                tokens.push_back(token);
            }
            if (tokens.size() < 5)
            {
                std::cerr << "SITE行格式错误: 需要至少5个字段，实际为 " << tokens.size() << std::endl;
                return 2;
            }
            Q = tokens[4].size();
            break;
        }
    }
    if (!found_site)
    {
        std::cerr << "未找到SITE行" << std::endl;
        return 2;
    }
    if (Q < 1)
    {
        std::cerr << "无效的M: " << Q << std::endl;
        return 3;
    }
    std::cerr << "M = " << Q << std::endl;

    // Step 2: 设置IDs
    qIDs.resize(Q);
    for (int i = 0; i < Q; i++)
    {
        qIDs[i] = "query_" + std::to_string(i);
    }

    // Step 3: 计算位点数 (N)
    int siteNum = 0;
    in.clear();
    in.seekg(0);
    while (std::getline(in, line))
    {
        if (line.rfind("SITE:", 0) == 0)
        {
            siteNum++;
        }
    }
    if (siteNum != N)
    {
        std::cerr << "panel跟query的位点数目不匹配 " << std::endl;
        return 4;
    }
    // Step 4: 初始化数据结构
    try
    {
        Z.resize(Q, std::vector<Syllable>(n));
        queryMultiInfo.resize(Q, std::vector<std::pair<unsigned int, uint8_t>>(n, std::make_pair(-1, 0)));

        querySyllableHavingMissing.resize(Q, vector<bool>(n));
    }
    catch (const std::bad_alloc& e)
    {
        std::cerr << "内存分配失败: " << e.what() << std::endl;
        return -1;
    }

    // Step 5: 处理SITE行
    in.clear();
    in.seekg(0);
    std::vector<Syllable> Z_(Q, 0);
    std::vector<std::vector<std::pair<uint8_t, uint8_t>>> syllableMultis(Q); // 临时存储多字符位点
    //Syllable filterTemp = 0;
    std::vector<Syllable> missingTemp(Q, 0); // 新增：query 的 missing 掩码
    Syllable one = 1;
    int K = 0, k = 0;
    while (std::getline(in, line))
    {
        if (line.rfind("SITE:", 0) != 0)
        {
            continue;
        }

        if (K >= N)
        {
            std::cerr << "SITE行数过多: K=" << K << ", 预期N=" << N << std::endl;
            return 10;
        }

        // 处理音节边界
        if (K % B == 0 && K != 0)
        {
            k = K / B - 1; // 上一个音节的索引
            if (k >= n)
            {
                std::cerr << "无效的k: " << k << ", n=" << n << std::endl;
                return 11;
            }
            for (int i = 0; i < Q; i++)
            {
                Z[i][k] = Z_[i]; // 保存压缩音节
                // 存储多字符位点
                if (!syllableMultis[i].empty())
                {
                    unsigned int start = queryMultiValues.size();
                    for (const auto& p : syllableMultis[i])
                    {
                        queryMultiValues.push_back(p);
                    }
                    queryMultiInfo[i][k] = std::make_pair(start, static_cast<uint8_t>(syllableMultis[i].size()));
                    syllableMultis[i].clear();
                }
                else
                {
                    queryMultiInfo[i][k] = std::make_pair(-1, 0);
                }

                if (missingTemp[i] != 0)
                {
                    queryMissingData[{i, k}] = missingTemp[i];
                    filter[k] = filter[k] | missingTemp[i];
                }
            }
            Z_.assign(Q, 0); // 重置X_
            missingTemp.assign(Q, 0);
        }

        std::stringstream ss(line);
        std::string token;
        std::getline(ss, token, '\t'); // Skip "SITE:"
        std::getline(ss, token, '\t'); // Skip index
        std::getline(ss, token, '\t'); // Skip physLoc
        std::getline(ss, token, '\t'); // Skip other column
        std::getline(ss, token, '\t'); // Get haplotype data

        if (token.size() != Q)
        {
            std::cerr << "单倍型数据长度不匹配: 预期 " << Q << ", 实际 " << token.size() << ", K=" << K << std::endl;
            return 6;
        }

        int index = 0;
        for (char c : token)
        {
            if (index >= Q)
            {
                std::cerr << "索引越界: index=" << index << ", M=" << M << ", K=" << K << std::endl;
                return 9;
            }
            if (c == '0')
            {
                Z_[index] = Z_[index] << 1;
                ++queryCount[0];
            }
            else if (c == '1')
            {
                Z_[index] = (Z_[index] << 1) | 1;
                ++queryCount[1];
            }
            else if (c > '1' && c <= '9')
            {
                Z_[index] = (Z_[index] << 1) | 1; // 多字符位点记为1
                ++queryCount[c-'0'];
                syllableMultis[index].push_back(
                    std::make_pair(static_cast<uint8_t>(K % B), static_cast<uint8_t>(c - '0')));
            }
            else if (c == '.')
            {
                Z_[index] = (Z_[index] << 1) | 1;
                ++queryCount[10];
                missingTemp[index] |= (one << (B - 1 - K % B)); // 记录 query missing
                querySyllableHavingMissing[index][K / B] = true;
            }
            else
            {
                std::cerr << "无效字符 '" << c << "' 在 K=" << K << ", index=" << index << std::endl;
                return 8;
            }
            //queryCount[c - '0'] += 1;
            index++;
        }
        if (index != Q)
        {
            std::cerr << "处理了 " << index << " 个单倍型，预期 " << Q << ", K=" << K << std::endl;
            return 9;
        }

        // 处理最后一个音节
        if (K == N - 1)
        {
            k = K / B; // 最后一个音节的索引
            int pad2 = n * B - N; // 填充位点数

            for (int i = 0; i < Q; i++)
            {
                Z_[i] <<= pad2; // 填充0
                Z[i][k] = Z_[i];
                if (!syllableMultis[i].empty())
                {
                    unsigned int start = queryMultiValues.size();
                    for (const auto& p : syllableMultis[i])
                    {
                        queryMultiValues.push_back(p);
                    }
                    queryMultiInfo[i][k] = std::make_pair(start, static_cast<uint8_t>(syllableMultis[i].size()));
                    syllableMultis[i].clear();
                }
                else
                {
                    queryMultiInfo[i][k] = std::make_pair(-1, 0);
                }
                if (querySyllableHavingMissing[i][k] == true)
                {
                    queryMissingData[{i, k}] = missingTemp[i];
                }
            }
        }
        K++;
    }

    if (K != N)
    {
        std::cerr << "处理了 " << K << " 个位点，预期 " << N << std::endl;
        return 10;
    }
    queryMultiValues.shrink_to_fit();
    end = clock();
    readQueryTime = ((double)(end - start)) / CLOCKS_PER_SEC;
    std::cerr << "readPanelTime = " << readPanelTime << " 秒, 多字符位点数 = " << panelMultiValues.size() << std::endl;
    return 0;
}

template <class Syllable>
int wmFSPBWT<Syllable>::makePanel() // fuzzy way :overall
{
    //get fuzzyX
    clock_t start, end;
    start = clock();
    if (FF * n % 32 != 0)
    {
        fuzzyX.resize(M, vector<uint32_t>(FF * n / 32 + 1));
    }
    else
    {
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
    for (int k = 0; k < n; k++)
    {
        if (count == 32 / FF)
        {
            ++index;
            count = 0;
        }
        for (int i = 0; i < M; i++)
        {
            int num = 0;
            if (B == 64)
            {
                num = __builtin_popcountll(X[i][k] | filter[k]);
            }
            else if (B == 128)
            {
                num = countSetBits128(X[i][k] | filter[k]);
            }
            uint32_t temp = num % T;
            fuzzyX[i][index] = (fuzzyX[i][index] << FF) | temp;
        }
        count++;
    }
    int pad1 = 32 - count * FF;
    for (int i = 0; i < M; i++)
    {
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
        for (int _ = 0; _ < T; _++)
        {
            p[_] = k + 1;
        }
        unit_index = k / (32 / FF);
        offset = 32 - FF - k % (32 / FF) * FF; //get fuzzy from fuzzyX

        for (int i = 0; i < M; i++) //each hapolotype
        {
            //update u[][][]
            for (int _ = 0; _ < T; _++)
            {
                u[k * (M * T) + i * T + _] = a_count[_];
                //u[k][i][_] = a_count[_];
                if (divergence[k][i] > p[_])
                {
                    p[_] = divergence[k][i];
                }
            }
            index1 = array[k][i];
            //uint32_t temp = fuzzyX[index][unit_index];

            temp = fuzzyX[index1][unit_index];
            if (F == 1)
            {
                fuzzy = (temp >> offset) & 0b1;
            }
            else if (F == 2)
            {
                fuzzy = (temp >> offset) & 0b11;
            }
            else if (F == 3)
            {
                fuzzy = (temp >> offset) & 0b111;
            }
            else if (F == 4)
            {
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
        for (int _ = 0; _ < T; _++)
        {
            for (int w = 0; w < a_count[_]; w++)
            {
                array[k + 1][m] = a[_][w];
                divergence[k + 1][m] = d[_][w];
                ++m;
            }
        }
        if (m != M)
        {
            return 4;
        }
        //put plus into u for w()
        for (int _ = 1; _ < T; _++)
        {
            plus = 0;
            for (int j = 0; j < _; j++)
            {
                plus += a_count[j];
            }
            for (int j = 0; j < M; j++)
            {
                u[k * (M * T) + j * T + _] += plus;
            }
        }
        for (int _ = 0; _ < T; _++)
        {
            a_count[_] = 0;
            d_count[_] = 0;
        }
    }
    end = clock();
    makeFuzzyPanelTime = ((double)(end - start)) / CLOCKS_PER_SEC;
    return 0;
}

template <class Syllable>
int wmFSPBWT<Syllable>::inPanelLongMatchQuery(int L, string inPanelOutput_file)
{
    clock_t start, end;
    start = clock();
    if (L < minSiteL)
    {
        return 1;
    }
    ofstream out(inPanelOutput_file);
    if (out.fail())
        return 2;

    int l = (L - (B - 1)) / B;
    int k;
    for (k = 0; k < n - 1; k++)
    {
        int unit_index = k / (32 / FF);
        int offset = 32 - FF - k % (32 / FF) * FF;

        bool m[T];
        for (int _ = 0; _ < T; _++)
        {
            m[_] = false;
        }
        int top = 0;
        bool report = false;
        for (int i = 0; i < M; i++)
        {
            if (divergence[k][i] > k - l)
            {
                for (int w = 0; w < T - 1; w++)
                {
                    for (int v = w + 1; v < T; v++)
                    {
                        if (m[w] == true && m[v] == true)
                        {
                            report = true;
                            break;
                        }
                    }
                }
                if (report == true)
                {
                    for (int i_a = top; i_a < i - 1; i_a++)
                    {
                        int maxDivergence = 0;
                        for (int i_b = i_a + 1; i_b < i; i_b++)
                        {
                            if (divergence[k][i_b] > maxDivergence)
                            {
                                maxDivergence = divergence[k][i_b];
                            }
                            uint32_t temp1, temp2, fuzzy1, fuzzy2;
                            int index_a = array[k][i_a], index_b =
                                    array[k][i_b];
                            temp1 = fuzzyX[index_a][unit_index];
                            temp2 = fuzzyX[index_b][unit_index];
                            if (F == 1)
                            {
                                fuzzy1 = (temp1 >> offset) & 0b1;
                                fuzzy2 = (temp2 >> offset) & 0b1;
                            }
                            else if (F == 2)
                            {
                                fuzzy1 = (temp1 >> offset) & 0b11;
                                fuzzy2 = (temp2 >> offset) & 0b11;
                            }
                            else if (F == 3)
                            {
                                fuzzy1 = (temp1 >> offset) & 0b111;
                                fuzzy2 = (temp2 >> offset) & 0b111;
                            }
                            else if (F == 4)
                            {
                                fuzzy1 = (temp1 >> offset) & 0b1111;
                                fuzzy2 = (temp2 >> offset) & 0b1111;
                            }

                            if (fuzzy1 != fuzzy2)
                            {
                                inPanelIdentification(L, maxDivergence - 1,
                                                      k, index_a, index_b, out);
                                ++inPanelAlternativeNum;
                            }
                        }
                    }
                    report = false;
                } //end if(report==true)

                top = i;
                for (int _ = 0; _ < T; _++)
                {
                    m[_] = false;
                }
                //report = false;
            } // end if(divergence[k][i]>k-l)
            //change m[]

            /*int fuzzy;
             fuzzy = catchFuzzy(array[k][i], k);*/

            uint32_t temp, fuzzy;
            temp = fuzzyX[array[k][i]][unit_index];
            if (F == 1)
            {
                fuzzy = (temp >> offset) & 0b1;
            }
            else if (F == 2)
            {
                fuzzy = (temp >> offset) & 0b11;
            }
            else if (F == 3)
            {
                fuzzy = (temp >> offset) & 0b111;
            }
            else if (F == 4)
            {
                fuzzy = (temp >> offset) & 0b1111;
            }

            m[fuzzy] = true;
        } //for i from 0 to M-1
        //cheak bottom block
        for (int w = 0; w < T - 1; w++)
        {
            for (int v = w + 1; v < T; v++)
            {
                if (m[w] == true && m[v] == true)
                {
                    report = true;
                }
            }
        }
        if (report == true)
        {
            for (int i_a = top; i_a < M - 1; i_a++)
            {
                int maxDivergence = 0;
                for (int i_b = i_a + 1; i_b < M; i_b++)
                {
                    if (divergence[k][i_b] > maxDivergence)
                    {
                        maxDivergence = divergence[k][i_b];
                    }
                    uint32_t temp1, temp2, fuzzy1, fuzzy2;
                    int index_a = array[k][i_a];
                    int index_b = array[k][i_b];
                    temp1 = fuzzyX[index_a][unit_index];
                    temp2 = fuzzyX[index_b][unit_index];
                    if (F == 1)
                    {
                        fuzzy1 = (temp1 >> offset) & 0b1;
                        fuzzy2 = (temp2 >> offset) & 0b1;
                    }
                    else if (F == 2)
                    {
                        fuzzy1 = (temp1 >> offset) & 0b11;
                        fuzzy2 = (temp2 >> offset) & 0b11;
                    }
                    else if (F == 3)
                    {
                        fuzzy1 = (temp1 >> offset) & 0b111;
                        fuzzy2 = (temp2 >> offset) & 0b111;
                    }
                    else if (F == 4)
                    {
                        fuzzy1 = (temp1 >> offset) & 0b1111;
                        fuzzy2 = (temp2 >> offset) & 0b1111;
                    }

                    if (fuzzy1 != fuzzy2)
                    {
                        inPanelIdentification(L, maxDivergence - 1, k,
                                              index_a, index_b, out);
                        ++inPanelAlternativeNum;
                    }
                }
            }
        }
    }
    int top = 0;
    for (int i = 0; i < M; i++)
    {
        if (divergence[k][i] > k - l + 1)
        {
            for (int i_a = top; i_a < i - 1; i_a++)
            {
                int maxDivergence = 0;
                for (int i_b = i_a + 1; i_b < i; i_b++)
                {
                    if (divergence[k][i_b] > maxDivergence)
                    {
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

                    if (F == 1)
                    {
                        fuzzy1 = (temp1 >> offset) & 0b1;
                        fuzzy2 = (temp2 >> offset) & 0b1;
                    }
                    else if (F == 2)
                    {
                        fuzzy1 = (temp1 >> offset) & 0b11;
                        fuzzy2 = (temp2 >> offset) & 0b11;
                    }
                    else if (F == 3)
                    {
                        fuzzy1 = (temp1 >> offset) & 0b111;
                        fuzzy2 = (temp2 >> offset) & 0b111;
                    }
                    else if (F == 4)
                    {
                        fuzzy1 = (temp1 >> offset) & 0b1111;
                        fuzzy2 = (temp2 >> offset) & 0b1111;
                    }

                    if (fuzzy1 == fuzzy2)
                    {
                        inPanelIdentification(L, maxDivergence - 1, n,
                                              index_a, index_b, out);
                        ++inPanelAlternativeNum;
                    }
                    else if (fuzzy1 != fuzzy2)
                    {
                        if (k - maxDivergence >= l)
                        {
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
    for (int i_a = top; i_a < M - 1; i_a++)
    {
        int maxDivergence = 0;
        for (int i_b = i_a + 1; i_b < M; i_b++)
        {
            int index_a = array[k][i_a];
            int index_b = array[k][i_b];
            if (divergence[k][i_b] > maxDivergence)
            {
                maxDivergence = divergence[k][i_b];
            }
            inPanelIdentification(L, maxDivergence - 1, n, index_a, index_b,
                                  out);
            ++inPanelAlternativeNum;
        }
    }

    end = clock();

    inPanelQueryTime = ((double)(end - start)) / CLOCKS_PER_SEC;

    inPanelFilteringTime = inPanelQueryTime - inPanelIdentificationTime;
    return 0;
    out.close();
    cout << "matches has been put into " << inPanelOutput_file << endl;
}

template <class Syllable>
int wmFSPBWT<Syllable>::outPanelLongMatchQuery(int L, string outPanelOutput_file)
{
    clock_t start, end;
    start = clock();

    if (FF * n % 32 != 0)
    {
        fuzzyZ.resize(Q, vector<uint32_t>(FF * n / 32 + 1));
    }
    else
    {
        fuzzyZ.resize(Q, vector<uint32_t>(FF * n / 32));
    }
    int index = 0, count = 0;
    for (int k = 0; k < n; k++)
    {
        if (count == 32 / FF)
        {
            ++index;
            count = 0;
        }
        for (int i = 0; i < Q; i++)
        {
            int num = 0;
            if (B == 64)
            {
                num = __builtin_popcountll((Z[i][k] | filter[k]));
            }
            else if (B == 128)
            {
                num = countSetBits128((Z[i][k] | filter[k]));
            }
            uint32_t temp = num % T;
            fuzzyZ[i][index] = (fuzzyZ[i][index] << FF) | temp;
        }
        count++;
    }
    int pad1 = 32 - count * FF;
    for (int i = 0; i < Q; i++)
    {
        fuzzyZ[i][fuzzyZ[0].size() - 1] =
            fuzzyZ[i][fuzzyZ[0].size() - 1] << pad1;
    }

    if (L < minSiteL)
    {
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
    for (int q = 0; q < Q; q++)
    {
        fill(dZ.begin(), dZ.end(), 0);
        fill(t.begin(), t.end(), 0);
        fill(Zdivergence.begin(), Zdivergence.end(), 0);
        fill(belowZdivergence.begin(), belowZdivergence.end(), 0);

        string qID = qIDs[q >> 1] + "-" + to_string(q & 1);
        t[0] = 0;

        //fake location
        for (int k = 0; k < n; k++)
        {
            uint32_t temp, fuzzy;
            int unit_index = k / (32 / FF);
            int offset = 32 - FF - k % (32 / FF) * FF;
            temp = fuzzyZ[q][unit_index];

            if (F == 1)
            {
                fuzzy = (temp >> offset) & 0b1;
            }
            else if (F == 2)
            {
                fuzzy = (temp >> offset) & 0b11;
            }
            else if (F == 3)
            {
                fuzzy = (temp >> offset) & 0b111;
            }
            else if (F == 4)
            {
                fuzzy = (temp >> offset) & 0b1111;
            }

            if (t[k] != M)
            {
                t[k + 1] =
                    u[k * (M * T) + t[k] * T + fuzzy];
                //                    u[k][t[k]][fuzzy];
            }
            else
            {
                // t[k] == M
                if (fuzzy < T - 1)
                {
                    t[k + 1] = u[k * (M * T) + 0 * T + fuzzy + 1];
                    //u[k][0][fuzzy + 1];
                }
                else if (fuzzy == T - 1)
                {
                    t[k + 1] = M;
                }
                else
                {
                    return 2;
                }
            }
        }

        Zdivergence[n + 1] = belowZdivergence[n + 1] = n;
        for (int k = n; k >= 0; --k)
        {
            Zdivergence[k] = std::min(Zdivergence[k + 1], k);
            belowZdivergence[k] = std::min(belowZdivergence[k + 1], k);
            if (t[k] != 0)
            {
                uint32_t panelTemp, panelFuzzy, queryTemp, queryFuzzy;
                int unit_index = (Zdivergence[k] - 1) / (32 / FF);
                int reminder = (Zdivergence[k] - 1) % (32 / FF);
                int offset = 32 - FF - reminder * FF;

                int index = array[k][t[k] - 1];
                //int index = array[Zdivergence[k] - 1][t[k] - 1];	//hapolotype on query
                panelTemp = fuzzyX[index][unit_index];
                queryTemp = fuzzyZ[q][unit_index];
                if (F == 1)
                {
                    panelFuzzy = (panelTemp >> offset) & 0b1;
                    queryFuzzy = (queryTemp >> offset) & 0b1;
                }
                else if (F == 2)
                {
                    panelFuzzy = (panelTemp >> offset) & 0b11;
                    queryFuzzy = (queryTemp >> offset) & 0b11;
                }
                else if (F == 3)
                {
                    panelFuzzy = (panelTemp >> offset) & 0b111;
                    queryFuzzy = (queryTemp >> offset) & 0b111;
                }
                else if (F == 4)
                {
                    panelFuzzy = (panelTemp >> offset) & 0b1111;
                    queryFuzzy = (queryTemp >> offset) & 0b1111;
                }
                //向前更新Zdivergence
                while (Zdivergence[k] > 0 && panelFuzzy == queryFuzzy)
                {
                    --Zdivergence[k];
                    if (reminder == 0)
                    {
                        if (unit_index == 0)
                        {
                            break;
                        }
                        else
                        {
                            reminder = 32 / FF - 1; //last bug is here ^^
                            --unit_index;
                        }
                    }
                    else
                    {
                        --reminder;
                    }
                    offset = 32 - FF - reminder * FF;
                    panelTemp = fuzzyX[index][unit_index];
                    queryTemp = fuzzyZ[q][unit_index];
                    if (F == 1)
                    {
                        panelFuzzy = (panelTemp >> offset) & 0b1;
                        queryFuzzy = (queryTemp >> offset) & 0b1;
                    }
                    else if (F == 2)
                    {
                        panelFuzzy = (panelTemp >> offset) & 0b11;
                        queryFuzzy = (queryTemp >> offset) & 0b11;
                    }
                    else if (F == 3)
                    {
                        panelFuzzy = (panelTemp >> offset) & 0b111;
                        queryFuzzy = (queryTemp >> offset) & 0b111;
                    }
                    else if (F == 4)
                    {
                        panelFuzzy = (panelTemp >> offset) & 0b1111;
                        queryFuzzy = (queryTemp >> offset) & 0b1111;
                    }
                }
            }
            else
            {
                //t[k]==0
                Zdivergence[k] = k;
            }
            if (t[k] < M)
            {
                uint32_t panelTemp, panelFuzzy, queryTemp, queryFuzzy;
                int unit_index = (belowZdivergence[k] - 1) / (32 / FF);
                int reminder = (belowZdivergence[k] - 1) % (32 / FF);
                int offset = 32 - FF - reminder * FF;
                int index = array[k][t[k]]; //hapolotype below query
                panelTemp = fuzzyX[index][unit_index];
                queryTemp = fuzzyZ[q][unit_index];
                if (F == 1)
                {
                    panelFuzzy = (panelTemp >> offset) & 0b1;
                    queryFuzzy = (queryTemp >> offset) & 0b1;
                }
                else if (F == 2)
                {
                    panelFuzzy = (panelTemp >> offset) & 0b11;
                    queryFuzzy = (queryTemp >> offset) & 0b11;
                }
                else if (F == 3)
                {
                    panelFuzzy = (panelTemp >> offset) & 0b111;
                    queryFuzzy = (queryTemp >> offset) & 0b111;
                }
                else if (F == 4)
                {
                    panelFuzzy = (panelTemp >> offset) & 0b1111;
                    queryFuzzy = (queryTemp >> offset) & 0b1111;
                }
                //向前更新belowZdivergence

                while (belowZdivergence[k] > 0 && panelFuzzy == queryFuzzy)
                {
                    belowZdivergence[k]--;
                    if (reminder == 0)
                    {
                        if (unit_index == 0)
                        {
                            break;
                        }
                        reminder = 32 / FF - 1;
                        --unit_index;
                    }
                    else
                    {
                        --reminder;
                    }

                    offset = 32 - FF - reminder * FF;
                    panelTemp = fuzzyX[index][unit_index];
                    queryTemp = fuzzyZ[q][unit_index];
                    if (F == 1)
                    {
                        panelFuzzy = (panelTemp >> offset) & 0b1;
                        queryFuzzy = (queryTemp >> offset) & 0b1;
                    }
                    else if (F == 2)
                    {
                        panelFuzzy = (panelTemp >> offset) & 0b11;
                        queryFuzzy = (queryTemp >> offset) & 0b11;
                    }
                    else if (F == 3)
                    {
                        panelFuzzy = (panelTemp >> offset) & 0b111;
                        queryFuzzy = (queryTemp >> offset) & 0b111;
                    }
                    else if (F == 4)
                    {
                        panelFuzzy = (panelTemp >> offset) & 0b1111;
                        queryFuzzy = (queryTemp >> offset) & 0b1111;
                    }
                }
            }
            else
            {
                belowZdivergence[k] = k;
            }
        }

        int f, g;
        f = g = t[0];
        vector<int> ftemp, gtemp;
        ftemp.resize(T);
        gtemp.resize(T);

        for (int k = 0; k < n; k++)
        {
            int unit_index = k / (32 / FF);
            int offset = 32 - FF - k % (32 / FF) * FF;
            uint32_t temp = fuzzyZ[q][unit_index];
            uint32_t fuzzyQ;
            if (F == 1)
            {
                fuzzyQ = (temp >> offset) & 0b1;
            }
            else if (F == 2)
            {
                fuzzyQ = (temp >> offset) & 0b11;
            }
            else if (F == 3)
            {
                fuzzyQ = (temp >> offset) & 0b111;
            }
            else if (F == 4)
            {
                fuzzyQ = (temp >> offset) & 0b1111;
            }
            if (g == M)
            {
                if (f == M)
                {
                    //update ftemp
                    for (int i = 0; i < T; i++)
                    {
                        if (fuzzyQ != i)
                        {
                            if (i != T - 1)
                            {
                                ftemp[i] = u[k * (M * T) + 0 * T + i + 1];
                                // u[k][0][i + 1];
                            }
                            else
                            {
                                ftemp[i] = M;
                            }
                        }
                    }

                    if (fuzzyQ != T - 1)
                    {
                        f = u[k * (M * T) + 0 * T + fuzzyQ + 1];
                        //u[k][0][fuzzyQ + 1];
                    }
                    else
                    {
                        f = M;
                    }
                }
                else //f!=M
                {
                    for (int i = 0; i < T; i++)
                    {
                        if (fuzzyQ != i)
                        {
                            ftemp[i] = u[k * (M * T) + f * T + i];
                            //u[k][f][i];
                        }
                    }
                    f = u[k * (M * T) + f * T + fuzzyQ];
                    //u[k][f][fuzzyQ];
                }
                //update gtemp and g
                for (int i = 0; i < T; i++)
                {
                    if (fuzzyQ != i)
                    {
                        if (i < T - 1)
                        {
                            gtemp[i] = u[k * (M * T) + 0 * T + i + 1];
                            //u[k][0][i + 1];
                        }
                        else
                        {
                            gtemp[i] = M;
                        }
                    }
                }
                if (fuzzyQ < T - 1)
                {
                    g = u[k * (M * T) + 0 * T + fuzzyQ + 1];
                    //u[k][0][fuzzyQ + 1];
                }
                else
                {
                    g = M;
                }
            }
            else //g!=M
            {
                for (int i = 0; i < T; i++)
                {
                    if (i != fuzzyQ)
                    {
                        ftemp[i] = u[k * (M * T) + f * T + i];
                        //u[k][f][i];
                        gtemp[i] = u[k * (M * T) + g * T + i];
                        //u[k][g][i];
                    }
                }
                f = u[k * (M * T) + f * T + fuzzyQ];
                //u[k][f][fuzzyQ];
                g = u[k * (M * T) + g * T + fuzzyQ];
                //u[k][g][fuzzyQ];
            }

            //output matches
            for (int i = 0; i < T; i++)
            {
                if (i != fuzzyQ)
                {
                    while (ftemp[i] != gtemp[i])
                    {
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

            if (f == g)
            {
                if (k + 1 - Zdivergence[k + 1] == l)
                {
                    --f;
                    dZ[array[k + 1][f]] = k + 1 - l;
                    //store divergence
                }

                //if (k + 1 - belowZdivergence[k + 1] == l)
                if (k + 1 - belowZdivergence[k + 1] == l)
                {
                    //store divergence
                    dZ[array[k + 1][g]] = k + 1 - l;
                    ++g;
                }
            }
            if (f != g)
            {
                while (divergence[k + 1][f] <= k + 1 - l)
                {
                    --f;
                    dZ[array[k + 1][f]] = k + 1 - l;
                }
                while (g < M && divergence[k + 1][g] <= k + 1 - l)
                {
                    dZ[array[k + 1][g]] = k + 1 - l;

                    ++g;
                }
            }
        }

        //mathces no ending at
        while (f != g)
        {
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

    outPanelQueryTime = ((double)(end - start)) / CLOCKS_PER_SEC;

    outPanelFilteringTime = outPanelQueryTime - outPanelIdentificationTime;
    cout << "matches has been put into " << outPanelOutput_file << endl;
    return 0;
}


template <class Syllable>
bool wmFSPBWT<Syllable>::inPanelSyllableEqual(int index_a, int index_b, int k)
{
    Syllable tempA = X[index_a][k];
    Syllable tempB = X[index_b][k];

    // 1. 快速检查 (注意：filter 可能会掩盖 0 vs 1 的差异，所以这一步只能用来快速排除，不能用来确信)
    if ((tempA | filter[k]) != (tempB | filter[k]))
    {
        return false;
    }

    // 2. 准备缺失数据掩码
    Syllable missA = 0, missB = 0;
    auto keyA = std::make_pair(index_a, k);
    auto keyB = std::make_pair(index_b, k);
    if (panelMissingData.count(keyA)) missA = panelMissingData[keyA];
    if (panelMissingData.count(keyB)) missB = panelMissingData[keyB];

    // 3. 准备多等位数据映射
    auto info_a = panelMultiInfo[index_a][k];
    auto info_b = panelMultiInfo[index_b][k];

    // 只有当存在多等位位点时才构建 map，稍微优化性能
    bool hasMultiA = info_a.second > 0;
    bool hasMultiB = info_b.second > 0;

    // 我们不需要构建完整的 map，直接在循环里查找或者构建一次皆可
    // 为了保持逻辑清晰，这里还是构建 map，但只在需要时查找
    std::unordered_map<uint8_t, uint8_t> valA, valB;
    if (hasMultiA)
    {
        for (unsigned int i = 0; i < info_a.second; ++i)
        {
            auto [pos, v] = panelMultiValues[info_a.first + i];
            valA[pos] = v;
        }
    }
    if (hasMultiB)
    {
        for (unsigned int i = 0; i < info_b.second; ++i)
        {
            auto [pos, v] = panelMultiValues[info_b.first + i];
            valB[pos] = v;
        }
    }

    // 4. 逐位精确比对
    for (int j = 0; j < B; ++j)
    {
        if (j >= N - k * B) break;

        int bit_shift = B - 1 - j;
        Syllable bitMask = (Syllable)1 << bit_shift;

        // A. 检查缺失值：任一方为 missing 则视为匹配 (continue)
        if ((missA & bitMask) || (missB & bitMask))
        {
            continue;
        }

        // B. 获取原始比特值 (0 或 1)
        int rawBitA = (tempA >> bit_shift) & 1;
        int rawBitB = (tempB >> bit_shift) & 1;

        // C. 基础比对：如果原始比特不同，肯定不匹配
        if (rawBitA != rawBitB)
        {
            return false; // 这里捕捉到了被 filter 掩盖的 0 vs 1 差异
        }

        // D. 进阶比对：如果都是 1，才需要检查是否是多等位基因 (如 1 vs 2)
        if (rawBitA == 1)
        {
            uint8_t realA = 1;
            uint8_t realB = 1;

            if (hasMultiA && valA.count(j)) realA = valA[j];
            if (hasMultiB && valB.count(j)) realB = valB[j];

            if (realA != realB)
            {
                return false; // 捕捉 1 vs 2, 2 vs 3 等
            }
        }
    }

    return true;
}

template <class Syllable>
bool wmFSPBWT<Syllable>::outPanelSyllableEqual(int index_panel, int index_query, int k)
{
    Syllable p = X[index_panel][k];

    Syllable q = Z[index_query][k];


    // 1. 准备缺失数据掩码

    Syllable missP = 0, missQ = 0;

    if (panelMissingData.count({index_panel, k})) missP = panelMissingData[{index_panel, k}];

    if (queryMissingData.count({index_query, k})) missQ = queryMissingData[{index_query, k}];


    // 2. 快速检查

    // 逻辑：构建局部 filter (两者的 missing 并集)。

    // 只要有一方是 missing，该位就被置为 1 (视为匹配)。

    // 只有当两者都不是 missing 且 bit 值不同时，这里才会返回 false。

    Syllable localFilter = missP | missQ;

    if ((p | localFilter) != (q | localFilter))
    {
        return false;
    }


    // 3. 准备多等位数据映射 (优化：只在有多等位位点时构建)

    auto info_p = panelMultiInfo[index_panel][k];

    auto info_q = queryMultiInfo[index_query][k];


    bool hasMultiP = info_p.second > 0;

    bool hasMultiQ = info_q.second > 0;


    std::unordered_map<uint8_t, uint8_t> valP, valQ;

    if (hasMultiP)
    {
        for (unsigned int i = 0; i < info_p.second; ++i)
        {
            valP[panelMultiValues[info_p.first + i].first] = panelMultiValues[info_p.first + i].second;
        }
    }

    if (hasMultiQ)
    {
        for (unsigned int i = 0; i < info_q.second; ++i)
        {
            valQ[queryMultiValues[info_q.first + i].first] = queryMultiValues[info_q.first + i].second;
        }
    }


    // 4. 逐位精确比对

    for (int j = 0; j < B; ++j)
    {
        if (k * B + j >= N) break;


        int bit_shift = B - 1 - j;

        Syllable bitMask = (Syllable)1 << bit_shift;


        // A. 检查缺失值：任一方为 missing 则视为匹配，跳过后续检查

        if ((missP & bitMask) || (missQ & bitMask))
        {
            continue;
        }


        // B. 获取原始比特值 (0 或 1)

        int rawBitP = (p >> bit_shift) & 1;

        int rawBitQ = (q >> bit_shift) & 1;


        // C. 基础比对：如果原始比特不同，肯定不匹配

        if (rawBitP != rawBitQ)
        {
            return false;
        }


        // D. 进阶比对：如果都是 1，才检查多等位基因 (如 1 vs 2)

        if (rawBitP == 1)
        {
            uint8_t realP = 1;

            uint8_t realQ = 1;


            if (hasMultiP && valP.count(j)) realP = valP[j];

            if (hasMultiQ && valQ.count(j)) realQ = valQ[j];


            if (realP != realQ)
            {
                return false;
            }
        }
    }


    return true;
}

template <typename InfoPair, typename ValueVector>
inline uint8_t getMultiAlleleValue(int bit_pos,
                                   const InfoPair& info,
                                   const ValueVector& multiValues)
{
    if (info.second == 0) return 1;
    for (unsigned int i = 0; i < info.second; ++i) {
        if (multiValues[info.first + i].first == bit_pos) {
            return multiValues[info.first + i].second;
        }
    }
    return 1;
}

template <class Syllable>
void wmFSPBWT<Syllable>::inPanelRefine(int L, int s_idx, int e_idx, int index_a, int index_b, std::ofstream& out)
{
    // =========================================================
    // 调试开关：锁定报错的样本对 (13 和 21)
    // =========================================================
    string idA = IDs[index_a];
    string idB = IDs[index_b];

    // ==========================================
    // 1. 计算 Start (向左)
    // ==========================================
    int start = 0;
    if (s_idx == -1)
    {
        start = 0;
    }
    else
    {
        auto info_a = panelMultiInfo[index_a][s_idx];
        auto info_b = panelMultiInfo[index_b][s_idx];
        bool hasMissing = false;
        Syllable missA = 0, missB = 0;
        if (panelMissingData.count({index_a, s_idx})) {missA = panelMissingData[{index_a, s_idx}];hasMissing=true;}
        if (panelMissingData.count({index_b, s_idx})) {missB = panelMissingData[{index_b, s_idx}];hasMissing=true;}
        // --- [Fast Path]: 没有多等位点 且 没有缺失值 ---
        if (!hasMissing && info_a.second == 0 && info_b.second == 0) {
            unsigned long tz;
            Syllable xor_val = X[index_a][s_idx] ^ X[index_b][s_idx];

            if (xor_val == 0) {
                tz = B; // 整个块都匹配
            } else {
                if (B == 64) {
                    tz = __builtin_ctzll(xor_val);
                } else if (B == 128) {
                    tz = ctz128_uint128(xor_val);
                }
            }
            start = (s_idx + 1) * B - tz;
        }
        else
        {
            int suffix_len = 0;
            for (int j = B - 1; j >= 0; --j)
            {
                int bit_shift = B - 1 - j;
                Syllable bitMask = ((Syllable)1) << bit_shift;

                bool is_missing = (missA & bitMask) || (missB & bitMask);
                int val_a = (X[index_a][s_idx] >> bit_shift) & 1;
                int val_b = (X[index_b][s_idx] >> bit_shift) & 1;

                bool is_match = true;
                if (!is_missing)
                {
                    if (val_a != val_b) is_match = false;
                    else if (val_a == 1)
                    {
                        // Multi-allelic check (simplified for debug view)
                        int ra = 1, rb = 1;
                        if (info_a.second > 0)
                        {
                            for (unsigned k = 0; k < info_a.second; ++k) if (panelMultiValues[info_a.first + k].first == j)
                            {
                                ra = panelMultiValues[info_a.first + k].second;
                                break;
                            }
                        }
                        if (info_b.second > 0)
                        {
                            for (unsigned k = 0; k < info_b.second; ++k) if (panelMultiValues[info_b.first + k].first == j)
                            {
                                rb = panelMultiValues[info_b.first + k].second;
                                break;
                            }
                        }
                        if (ra != rb) is_match = false;
                    }
                }

                if (!is_match) break;
                suffix_len++;
            }
            start = (s_idx + 1) * B - suffix_len;
        }
    }


    int end = 0;

    if (e_idx == n) {
        end = N; // 已经匹配到基因组末尾
    } else {
        auto& info_a = panelMultiInfo[index_a][e_idx];
        auto& info_b = panelMultiInfo[index_b][e_idx];

        bool hasMissing = false;
        if (!panelMissingData.empty()) {
            if (panelMissingData.count({index_a, e_idx}) || panelMissingData.count({index_b, e_idx})) hasMissing = true;
        }

        // [Fast Path]
        if (!hasMissing && info_a.second == 0 && info_b.second == 0) {
            Syllable diff = X[index_a][e_idx] ^ X[index_b][e_idx];
            unsigned long lz = 0;
            // 理论上 diff 不会为 0，因为 Identification 在这里停下了
            // 但如果为 0，说明这个块也是匹配的（边界情况），lz = B
            if (diff == 0) lz = B;
            else if (B == 64) lz = __builtin_clzll(diff);
            else if (B == 128) lz = clz128_uint128(diff);

            end = e_idx * B + lz;
        }
        // [Slow Path]
        else {
            Syllable missA = 0, missB = 0;
            if (hasMissing) {
                if (panelMissingData.count({index_a, e_idx})) missA = panelMissingData[{index_a, e_idx}];
                if (panelMissingData.count({index_b, e_idx})) missB = panelMissingData[{index_b, e_idx}];
            }

            int match_len = 0;
            for (int j = 0; j < B; ++j) {
                int global_pos = e_idx * B + j;
                if (global_pos >= N) break;

                int bit_shift = B - 1 - j;
                Syllable bitMask = ((Syllable)1) << bit_shift;

                if ((missA & bitMask) || (missB & bitMask)) { match_len++; continue; }

                int val_a = (X[index_a][e_idx] >> bit_shift) & 1;
                int val_b = (X[index_b][e_idx] >> bit_shift) & 1;

                if (val_a != val_b) break; // Mismatch found
                if (val_a == 1) {
                    uint8_t ra = getMultiAlleleValue(j, info_a, panelMultiValues);
                    uint8_t rb = getMultiAlleleValue(j, info_b, panelMultiValues);
                    if (ra != rb) break;
                }
                match_len++;
            }
            end = e_idx * B + match_len;
        }
    }




    if (end - start >= L)
    {
        out << IDs[index_a] << '\t' << IDs[index_b] << '\t' << start << '\t' << end - 1 << '\n';
        ++inPanelMatchNum;
        matchLen += (end - start);
    }
}

template <class Syllable>
void wmFSPBWT<Syllable>::outPanelRefine(int L, int s_idx, int e_idx, int index_panel, int index_query, std::ofstream& out)
{
    // ==========================================
    // 1. 计算 Start
    // ==========================================
    int start = 0;
    if (s_idx == -1) {
        start = 0;
    } else {
        auto& info_p = panelMultiInfo[index_panel][s_idx];
        auto& info_q = queryMultiInfo[index_query][s_idx];

        bool hasMissing = false;
        if (!panelMissingData.empty() || !queryMissingData.empty()) {
            if (panelMissingData.count({index_panel, s_idx}) || queryMissingData.count({index_query, s_idx})) hasMissing = true;
        }

        // [Fast Path]
        if (!hasMissing && info_p.second == 0 && info_q.second == 0) {
            Syllable diff = X[index_panel][s_idx] ^ Z[index_query][s_idx];
            unsigned long tz = 0;
            if (diff == 0) tz = B;
            else if (B == 64) tz = __builtin_ctzll(diff);
            else if (B == 128) tz = ctz128_uint128(diff);

            start = (s_idx + 1) * B - tz;
        }
        // [Slow Path]
        else {
            int suffix_len = 0;
            Syllable missP = 0, missQ = 0;
            if (hasMissing) {
                if (panelMissingData.count({index_panel, s_idx})) missP = panelMissingData[{index_panel, s_idx}];
                if (queryMissingData.count({index_query, s_idx})) missQ = queryMissingData[{index_query, s_idx}];
            }
            for (int j = B - 1; j >= 0; --j) {
                int bit_shift = B - 1 - j;
                Syllable bitMask = ((Syllable)1) << bit_shift;

                if ((missP & bitMask) || (missQ & bitMask)) { suffix_len++; continue; }

                int val_p = (X[index_panel][s_idx] >> bit_shift) & 1;
                int val_q = (Z[index_query][s_idx] >> bit_shift) & 1;

                if (val_p != val_q) break;
                if (val_p == 1) {
                    uint8_t rp = getMultiAlleleValue(j, info_p, panelMultiValues);
                    uint8_t rq = getMultiAlleleValue(j, info_q, queryMultiValues);
                    if (rp != rq) break;
                }
                suffix_len++;
            }
            start = (s_idx + 1) * B - suffix_len;
        }
    }

    // ==========================================
    // 2. 计算 End (只看 e_idx)
    // ==========================================
    int end = 0;
    if (e_idx == n) {
        end = N;
    } else {
        auto& info_p = panelMultiInfo[index_panel][e_idx];
        auto& info_q = queryMultiInfo[index_query][e_idx];

        bool hasMissing = false;
        if (!panelMissingData.empty() || !queryMissingData.empty()) {
            if (panelMissingData.count({index_panel, e_idx}) || queryMissingData.count({index_query, e_idx})) hasMissing = true;
        }

        // [Fast Path]
        if (!hasMissing && info_p.second == 0 && info_q.second == 0) {
            Syllable diff = X[index_panel][e_idx] ^ Z[index_query][e_idx];
            unsigned long lz = 0;
            if (diff == 0) lz = B; // 理论不应发生，除非 logic bug 或 e_idx 越界
            else if (B == 64) lz = __builtin_clzll(diff);
            else if (B == 128) lz = clz128_uint128(diff);

            end = e_idx * B + lz;
        }
        // [Slow Path]
        else {
            Syllable missP = 0, missQ = 0;
            if (hasMissing) {
                if (panelMissingData.count({index_panel, e_idx})) missP = panelMissingData[{index_panel, e_idx}];
                if (queryMissingData.count({index_query, e_idx})) missQ = queryMissingData[{index_query, e_idx}];
            }
            int match_len = 0;
            for (int j = 0; j < B; ++j) {
                int global_pos = e_idx * B + j;
                if (global_pos >= N) break;

                int bit_shift = B - 1 - j;
                Syllable bitMask = ((Syllable)1) << bit_shift;

                if ((missP & bitMask) || (missQ & bitMask)) { match_len++; continue; }

                int val_p = (X[index_panel][e_idx] >> bit_shift) & 1;
                int val_q = (Z[index_query][e_idx] >> bit_shift) & 1;

                if (val_p != val_q) break;
                if (val_p == 1) {
                    uint8_t rp = getMultiAlleleValue(j, info_p, panelMultiValues);
                    uint8_t rq = getMultiAlleleValue(j, info_q, queryMultiValues);
                    if (rp != rq) break;
                }
                match_len++;
            }
            end = e_idx * B + match_len;
        }
    }

    if (end - start >= L) {
        out << IDs[index_panel] << '\t' << qIDs[index_query] << '\t' << start << '\t' << end - 1 << '\n';
        ++outPanelMatchNum;
        matchLen += (end - start);
    }
}
template <class Syllable>
void wmFSPBWT<Syllable>::inPanelIdentification(int L, int s_idx, int e_idx, int index_a,
                                               int index_b, std::ofstream& out)
{
    // // ==================== 调试信息：谁、在什么时候、用什么范围调用了我 ====================
    // std::cerr << "\n>>> inPanelIdentification 被调用！\n";
    // std::cerr << "    Hap" << index_a << " vs Hap" << index_b
    //           << "   L=" << L
    //           << "   s_idx=" << s_idx
    //           << "   e_idx=" << e_idx
    //           << "   (范围长度 = " << (e_idx - s_idx) << " 个 syllable)\n";
    // // =====================================================================================
    alternativeSyllableNum += (e_idx - s_idx + 1);
    clock_t start, end;
    start = clock();
    int l = (L - 2 * (B - 1)) / B;
    int head = s_idx + 1, tail;
    while (head < e_idx)
    {
        tail = head;
        while (tail < e_idx && inPanelSyllableEqual(index_a, index_b, tail))
        {
            tail++;
        }
        if (tail - head >= l)
        {
            inPanelRefine(L, head - 1, tail, index_a, index_b, out);
        }
        head = tail + 1;
        while (head < e_idx && !inPanelSyllableEqual(index_a, index_b, head))
        {
            head++;
        }
    }
    end = clock();
    inPanelIdentificationTime += ((double)(end - start)) / CLOCKS_PER_SEC;
}

template <class Syllable>
void wmFSPBWT<Syllable>::outPanelIdentification(int L, int s_idx, int e_idx, int index_panel,
                                                int index_query, std::ofstream& out)
{
    //
    // // ==================== 调试信息：谁、在什么时候、用什么范围调用了我 ====================
    // std::cerr << "\n>>> outPanelIdentification 被调用！\n";
    // std::cerr << " PanelHap" << index_panel << " vs QueryHap" << index_query
    //           << " L=" << L
    //           << " s_idx=" << s_idx
    //           << " e_idx=" << e_idx
    //           << " (范围长度 = " << (e_idx - s_idx) << " 个 syllable)\n";
    // // =====================================================================================


    alternativeSyllableNum += (e_idx - s_idx + 1);
    clock_t start, end;
    start = clock();
    int l = (L - 2 * (B - 1)) / B;
    int head = s_idx + 1, tail;
    while (head < e_idx)
    {
        tail = head;
        while (tail < e_idx && outPanelSyllableEqual(index_panel, index_query, tail))
        {
            tail++;
        }
        if (tail - head >= l)
        {
            outPanelRefine(L, head - 1, tail, index_panel, index_query, out);
        }
        head = tail + 1;
        while (head < e_idx && !outPanelSyllableEqual(index_panel, index_query, head))
        {
            head++;
        }
    }
    end = clock();
    outPanelIdentificationTime += ((double)(end - start)) / CLOCKS_PER_SEC;
}

template <class Syllable>
void wmFSPBWT<Syllable>::outputInformationToFile(const std::string& fileName, string mode)
{
    std::ofstream outputFile(fileName);

    if (outputFile.is_open())
    {
        if (mode == "in")
        {
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
            outputFile << "AlternativeLen: " << alternativeSyllableNum * B
                << std::endl;
            outputFile << "MatchLen: " << matchLen
                << std::endl;
            // 新增：输出面板位点计数信息
            outputFile << "Panel Site Counts:" << std::endl;
            for (int i = 0; i < 10; i++)
            {
                outputFile << "Site " << i << ": " << panelCount[i] << std::endl;
            }
            outputFile << "PanelMissingSite "  << panelCount[10] << std::endl;
            std::cout << "Information has been written to " << fileName
                << std::endl;
        }
        else if (mode == "out")
        {
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
            outputFile << "AlternativeLen: " << alternativeSyllableNum * B
                << std::endl;
            outputFile << "MatchLen: " << matchLen
                << std::endl;
            // 新增：输出面板和查询位点计数信息
            outputFile << "Panel Site Counts:" << std::endl;
            for (int i = 0; i < 10; i++)
            {
                outputFile << "Site " << i << ": " << panelCount[i] << std::endl;
            }
            outputFile << "PanelMissingSite "  << panelCount[10] << std::endl;
            outputFile << std::endl;

            outputFile << "Query Site Counts:" << std::endl;
            for (int i = 0; i < 10; i++)
            {
                outputFile << "Site " << i << ": " << queryCount[i] << std::endl;
            }
            outputFile << "QueryMissingSite "  << queryCount[10] << std::endl;
            outputFile << std::endl;
            outputFile.close();
            std::cout << "Information has been written to " << fileName
                << std::endl;
        }
    }
    else
    {
        std::cerr << "Unable to open file " << fileName << " for writing."
            << std::endl;
    }
}
template <class Syllable>
void wmFSPBWT<Syllable>::save(const string& index_file)
{
    std::ofstream out(index_file, std::ios::binary);
    if (!out.is_open()) {
        throw std::runtime_error("无法打开索引文件进行写入: " + index_file);
    }
    uint64_t magic = 0x574D4650'42425400ULL;  // "WMFSPBWT\0"
    out.write(reinterpret_cast<const char*>(&magic), sizeof(magic));

    uint32_t sz_size = sizeof(size_t);
    uint32_t syl_size = sizeof(Syllable);
    out.write(reinterpret_cast<const char*>(&sz_size), sizeof(sz_size));
    out.write(reinterpret_cast<const char*>(&syl_size), sizeof(syl_size));
    // 1. 保存基本标量成员
    out.write(reinterpret_cast<const char*>(&B), sizeof(B));
    out.write(reinterpret_cast<const char*>(&F), sizeof(F));
    out.write(reinterpret_cast<const char*>(&T), sizeof(T));
    out.write(reinterpret_cast<const char*>(&M), sizeof(M));
    out.write(reinterpret_cast<const char*>(&N), sizeof(N));
    out.write(reinterpret_cast<const char*>(&n), sizeof(n));
    out.write(reinterpret_cast<const char*>(&minSiteL), sizeof(minSiteL));

    // 2. 保存 IDs (vector<string>)
    size_t ids_size = IDs.size();
    out.write(reinterpret_cast<const char*>(&ids_size), sizeof(ids_size));
    for (const auto& id : IDs) {
        size_t len = id.length();
        out.write(reinterpret_cast<const char*>(&len), sizeof(len));
        out.write(id.c_str(), len);
    }

    // 3. 保存 physLocs (vector<int>)
    saveVector(out, physLocs);

    // 4. 保存 X (vector<vector<Syllable>>)
    saveNestedVector(out, X);

    // 5. 保存 fuzzyX (vector<vector<uint32_t>>)
    saveNestedVector(out, fuzzyX);

    // 6. 保存 array (vector<vector<int>>)
    saveNestedVector(out, array);

    // 7. 保存 divergence (vector<vector<int>>)
    saveNestedVector(out, divergence);

    // 8. 保存 filter (vector<Syllable>)
    saveVector(out, filter);

    // 9. 保存 panelSyllableHavingMissing (vector<vector<bool>>)
    {
        size_t outer = panelSyllableHavingMissing.size();
        out.write(reinterpret_cast<const char*>(&outer), sizeof(outer));
        for (const auto& inner : panelSyllableHavingMissing) {
            size_t inner_sz = inner.size();
            out.write(reinterpret_cast<const char*>(&inner_sz), sizeof(inner_sz));
            for (bool b : inner) {
                out.write(reinterpret_cast<const char*>(&b), sizeof(bool));
            }
        }
    }

    // 10. 保存 panelMultiInfo (vector<vector<pair<unsigned int, uint8_t>>>)
    saveNestedVector(out, panelMultiInfo);

    // 11. 保存 panelMultiValues (vector<pair<uint8_t, uint8_t>>)
    saveVector(out, panelMultiValues);

    // 12. 保存 panelMissingData (unordered_map<pair<int,int>, Syllable>)
    size_t missing_size = panelMissingData.size();
    out.write(reinterpret_cast<const char*>(&missing_size), sizeof(missing_size));
    for (const auto& [key, value] : panelMissingData) {
        out.write(reinterpret_cast<const char*>(&key.first), sizeof(key.first));
        out.write(reinterpret_cast<const char*>(&key.second), sizeof(key.second));
        out.write(reinterpret_cast<const char*>(&value), sizeof(value));
    }

    // 13. 保存 u (int* 动态数组)
    size_t u_total_size = static_cast<size_t>(n) * M * T;
    out.write(reinterpret_cast<const char*>(&u_total_size), sizeof(u_total_size));
    if (u_total_size > 0 && u != nullptr) {
        out.write(reinterpret_cast<const char*>(u), u_total_size * sizeof(int));
    }

    // 14. 保存 panelCount (固定大小数组 u_long panelCount[11])
    out.write(reinterpret_cast<const char*>(panelCount), sizeof(panelCount));

    // 注意：这里没有保存 query 相关的成员（Z, fuzzyZ, queryMulti* 等），因为 save 只用于 panel 索引

    out.close();

    if (!out.good()) {
        throw std::runtime_error("写入索引文件时发生错误: " + index_file);
    }

    std::cout << "索引已成功保存到: " << index_file << std::endl;
}
template <class Syllable>
void wmFSPBWT<Syllable>::load(const std::string& index_file) {
    std::ifstream in(index_file, std::ios::binary);
    if (!in.is_open()) {
        throw std::runtime_error("无法打开索引文件进行读取: " + index_file);
    }


    // 读取 magic number 和版本信息（建议在 save 中也添加对应写入）
    uint64_t magic;
    in.read(reinterpret_cast<char*>(&magic), sizeof(magic));
    if (!in.good()) throw std::runtime_error("读取 magic number 失败");
    if (magic != 0x574D4650'42425400ULL) {  // "WMFSPBWT\0"
        throw std::runtime_error("无效的索引文件格式（magic number 不匹配）");
    }
    std::cout << "[DEBUG] Magic number 校验通过" << std::endl;

    // 标量字段
    in.read(reinterpret_cast<char*>(&B), sizeof(B));
    if (!in.good()) throw std::runtime_error("读取 B 失败");


    in.read(reinterpret_cast<char*>(&F), sizeof(F));


    in.read(reinterpret_cast<char*>(&T), sizeof(T));


    in.read(reinterpret_cast<char*>(&M), sizeof(M));


    in.read(reinterpret_cast<char*>(&N), sizeof(N));


    in.read(reinterpret_cast<char*>(&n), sizeof(n));


    in.read(reinterpret_cast<char*>(&minSiteL), sizeof(minSiteL));


    // IDs (vector<string>)
    size_t ids_sz;
    in.read(reinterpret_cast<char*>(&ids_sz), sizeof(ids_sz));
    if (!in.good()) throw std::runtime_error("读取 IDs size 失败");


    IDs.resize(ids_sz);
    for (size_t i = 0; i < ids_sz; ++i) {
        size_t len;
        in.read(reinterpret_cast<char*>(&len), sizeof(len));
        if (!in.good()) throw std::runtime_error("读取 ID[" + std::to_string(i) + "] length 失败");
        IDs[i].resize(len);
        if (len > 0) {
            in.read(&IDs[i][0], len);
            if (!in.good()) throw std::runtime_error("读取 ID[" + std::to_string(i) + "] data 失败");
        }
    }

    // physLocs (vector<int>)
    size_t phys_sz;
    in.read(reinterpret_cast<char*>(&phys_sz), sizeof(phys_sz));
    if (!in.good()) throw std::runtime_error("读取 physLocs size 失败");

    physLocs.resize(phys_sz);
    if (phys_sz > 0) {
        in.read(reinterpret_cast<char*>(physLocs.data()), phys_sz * sizeof(int));
        if (!in.good()) throw std::runtime_error("读取 physLocs data 失败");
    }

    // X (vector<vector<Syllable>>)
    size_t x_outer;
    in.read(reinterpret_cast<char*>(&x_outer), sizeof(x_outer));
    if (!in.good()) throw std::runtime_error("读取 X outer size 失败");

    X.resize(x_outer);
    for (auto& row : X) {
        loadVector(in, row);
    }

    // fuzzyX
    size_t fuzzy_outer;
    in.read(reinterpret_cast<char*>(&fuzzy_outer), sizeof(fuzzy_outer));
    if (!in.good()) throw std::runtime_error("读取 fuzzyX outer size 失败");

    fuzzyX.resize(fuzzy_outer);
    for (auto& row : fuzzyX) {
        loadVector(in, row);
    }

    // array
    size_t array_outer;
    in.read(reinterpret_cast<char*>(&array_outer), sizeof(array_outer));
    if (!in.good()) throw std::runtime_error("读取 array outer size 失败");
    array.resize(array_outer);
    for (auto& row : array) {
        loadVector(in, row);
    }

    // divergence
    size_t div_outer;
    in.read(reinterpret_cast<char*>(&div_outer), sizeof(div_outer));
    if (!in.good()) throw std::runtime_error("读取 divergence outer size 失败");
    divergence.resize(div_outer);
    for (auto& row : divergence) {
        loadVector(in, row);
    }

    // filter
    loadVector(in, filter);

    // panelSyllableHavingMissing (vector<vector<bool>>)
    size_t outer_sz;
    in.read(reinterpret_cast<char*>(&outer_sz), sizeof(outer_sz));
    if (!in.good()) throw std::runtime_error("读取 panelSyllableHavingMissing outer size 失败");
    panelSyllableHavingMissing.resize(outer_sz);
    for (auto& inner : panelSyllableHavingMissing) {
        size_t inner_sz;
        in.read(reinterpret_cast<char*>(&inner_sz), sizeof(inner_sz));
        if (!in.good()) throw std::runtime_error("读取 inner size 失败");
        inner.resize(inner_sz);
        for (size_t i = 0; i < inner_sz; ++i) {
            bool value;
            in.read(reinterpret_cast<char*>(&value), sizeof(bool));
            if (!in.good()) throw std::runtime_error("读取 bool 值失败");
            inner[i] = value;
        }
    }

    // panelMultiInfo
    size_t multiinfo_outer;
    in.read(reinterpret_cast<char*>(&multiinfo_outer), sizeof(multiinfo_outer));
    panelMultiInfo.resize(multiinfo_outer);
    for (auto& row : panelMultiInfo) {
        loadVector(in, row);
    }

    // panelMultiValues
    loadVector(in, panelMultiValues);

    // panelMissingData
    size_t missing_sz;
    in.read(reinterpret_cast<char*>(&missing_sz), sizeof(missing_sz));
    if (!in.good()) throw std::runtime_error("读取 panelMissingData size 失败");
    panelMissingData.clear();
    for (size_t i = 0; i < missing_sz; ++i) {
        std::pair<int, int> key;
        Syllable val;
        in.read(reinterpret_cast<char*>(&key.first), sizeof(key.first));
        in.read(reinterpret_cast<char*>(&key.second), sizeof(key.second));
        in.read(reinterpret_cast<char*>(&val), sizeof(val));
        if (!in.good()) throw std::runtime_error("读取 panelMissingData 条目失败");
        panelMissingData[key] = val;
    }

    // u (动态数组)
    size_t u_sz;
    in.read(reinterpret_cast<char*>(&u_sz), sizeof(u_sz));
    if (!in.good()) throw std::runtime_error("读取 u size 失败");


    if (u) delete[] u;
    u = new int[u_sz];
    if (u_sz > 0) {
        in.read(reinterpret_cast<char*>(u), u_sz * sizeof(int));
        if (!in.good()) throw std::runtime_error("读取 u data 失败");
    }

    // panelCount
    in.read(reinterpret_cast<char*>(panelCount), sizeof(panelCount));
    if (!in.good()) throw std::runtime_error("读取 panelCount 失败");

    in.close();
    if (!in.good()) {
        throw std::runtime_error("读取索引文件时发生错误（文件可能被截断）");
    }
    std::cout << "索引加载成功: " << index_file << std::endl;
}