
#include "cuda_nblist.h"
#include "doctest/doctest.h"
#include <fstream>
#include <vector>
#include <string>
#include <sstream>

namespace dpnblist {
 
    void read_lmp(std::string filename, std::vector<std::vector<float>> &xyz, Vec3<float> &box_len){
        
        std::ifstream file(filename);

        if (file.is_open()) {
            std::string line;

            // read atom numbers
            for(int i = 0; i < 3; ++i){
                std::getline(file, line);
            }
            std::istringstream atomLineStream(line);
            int numAtoms;
            atomLineStream >> numAtoms;
            xyz.resize(numAtoms,std::vector<float>(3,0.0f));
            
            // read box size
            while (std::getline(file, line) && line.find("xlo xhi") == std::string::npos);
            
            float lo, hi;
            for (int i = 0; i < 3; ++i) {
                std::istringstream boxLineStream(line);
                std::getline(file, line);
                boxLineStream >> lo >> hi;
                box_len[i] = hi - lo;
            }

            // jump exter lines
            while (std::getline(file, line) && line.find("Atoms") == std::string::npos);
            std::getline(file, line);
            // while (std::getline(file, line) && !line.empty());

            // read atom coordinates
            int id, type;
            float x, y, z;
            for (int i = 0; i < numAtoms; ++i) {
                std::getline(file, line);
                std::istringstream iss(line);
                iss >> id >> type >> x >> y >> z;
                xyz[i][0] = x;
                xyz[i][1] = y;
                xyz[i][2] = z;
            }

            file.close();
        } else {
            std::cerr << "Unable to open file: " << filename << std::endl;
        }
    }

    void read_ref(std::string filename, std::vector<std::vector<size_t>> &ref_listArray){
            std::ifstream file(filename);
            if (!file.is_open()) {
                std::cout << "无法打开文件." << std::endl;
            }
            
            std::string line;
            int n = 0;
            while (std::getline(file, line)) {
                std::vector<std::string> tokens;
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
            file.close();
    }

    void nblist_test(std::string filename, std::string reffilename, float cutoff, float skin = 0.0){
        std::vector<std::vector<float>> xyz;
        Vec3<float> box_len;
        read_lmp(filename, xyz, box_len);
        
        Box* box = new Box(box_len);
        CudaCellList nblist(box, cutoff, skin);
        nblist.build(xyz);
        nblist.update(xyz);
        std::vector<std::vector<size_t>> listArray = nblist.get_listArray();
        // Open a file ref file
        std::vector<std::vector<size_t>> ref_listArray(listArray.size());
        read_ref(reffilename, ref_listArray);
        CHECK_EQ(listArray.size(), ref_listArray.size());
        for (int i = 0; i < listArray.size(); ++i) {
            CHECK_EQ(listArray[i].size(), ref_listArray[i].size());
            if (listArray[i].size() != ref_listArray[i].size()) {
                std::cout << filename << " line wrong: " << i+1 << std::endl;
            }
        }
    }

    TEST_SUITE("Test Ref Cell List")
    {

        TEST_CASE("CudaCellList"){
        // test1
            nblist_test("heterogeneous.lmp", "ref_nblist_heterogeneous.txt", 3.0, 0.0);

        // test2
            // nblist_test("homogeneous_sparse.lmp", "sort_output_homogeneous_sparse.txt", 3.0, 0.0);

        // test3
            nblist_test("homogeneous_dense.lmp", "ref_nblist_homogeneous_dense.txt", 3.0, 0.0);

        // test4
            nblist_test("Nb_water.lmp", "ref_nblist_Nb_water.txt", 3.6, 0.0);

        // test5
            nblist_test("protein.lmp", "ref_nblist_protein.txt", 3.6, 0.0);

        }

    }

}
