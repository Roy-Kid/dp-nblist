#include <chemfiles.hpp>

#include "refCellList.h"

namespace dpnblist {

    TEST_SUITE("Test Ref Cell List")
    {

        TEST_CASE("NeighborList"){
            Box box({80.1, 80.05, 80.13});
            NeighborList nblist(&box, 2.6);
            std::vector<Vec3<double>> xyz;
            xyz.resize(50000);

            std::ifstream file1("50000.pdb");
            if (!file1.is_open()) {
                std::cout << "无法打开文件." << std::endl;
            }
            std::string line;
            // 跳过前5行
            for(int i=0;i<5;++i){
                std::getline(file1, line);
            }
            // 记录原子序号
            int n = 0;
            while (std::getline(file1, line)) {
                std::vector<std::string> tokens;
                // 使用字符串流将每行内容拆分
                std::stringstream ss(line);
                std::string token;
                while (ss >> token) {
                    tokens.push_back(token);
                }
                int len = tokens.size();
                if(len>1){
                    xyz[n][0] = stof(tokens[len-5]);
                    xyz[n][1] = stof(tokens[len-4]);
                    xyz[n][2] = stof(tokens[len-3]);
                }
                n++;
            }
            file1.close();
            // 测试build
            nblist.build(xyz);
            std::vector<std::vector<size_t>> listArray = nblist.get_listArray();
            // Open a file for writing
            std::ofstream outFile("output.txt");
            if (outFile.is_open()) {
                for (const auto& row : listArray) {
                    for (size_t num : row) {
                        outFile << num << '\t';
                    }
                    outFile << '\n';
                }

                outFile.close();
                std::cout << "File writing complete." << std::endl;
            } else {
                std::cerr << "Unable to open the file for writing." << std::endl;
            }

            std::vector<std::vector<size_t>> ref_listArray(50000);
            std::ifstream file2("out_nbl.txt");
            if (!file2.is_open()) {
                std::cout << "无法打开文件." << std::endl;
            }
            // 记录原子序号
            n = 0;
            while (std::getline(file2, line)) {
                std::vector<std::string> tokens;
                // 使用字符串流将每行内容拆分
                std::stringstream ss(line);
                std::string token;
                while (ss >> token) {
                    tokens.push_back(token);
                }
                int len = tokens.size();
                if(len>1){
                    for(int i = 1; i < len; ++i){
                        ref_listArray[n].push_back(stof(tokens[i]));
                    }
                    n++;
                }
            }
            file2.close();

            CHECK(listArray == ref_listArray);
        }

    }

}
