/*
Copyright(c) 2013, Ilya Vorobyev und Vasiliy Usatyuk
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met :
*Redistributions of source code must retain the above copyright
notice, this list of conditions and the following disclaimer.
* Redistributions in binary form must reproduce the above copyright
notice, this list of conditions and the following disclaimer in the
documentation and / or other materials provided with the distribution.
* Neither the name of the <organization> nor the
names of its contributors may be used to endorse or promote products
derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED.IN NO EVENT SHALL <COPYRIGHT HOLDER> BE LIABLE FOR ANY
DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/
#include".\myLib\TBlocalOptimization.h"
#include".\myLib\irregularLDPC.h"



void test() {
    time_t start = time(NULL);
    int variableNodes = 7;
    int checkNodes = 3;
    vector<vector<int> > protograph(3, vector<int>(7, 1));
    ll iterationCount = 0;
    ll circulant = 230;
    ll seed = 85;
    srand(seed);
    int count = 0;
    while (true) {
        ++count;
        vector<vector<vector<int> > > a(checkNodes, vector<vector<int> >(variableNodes));
        for (int r = 0; r < checkNodes; ++r) {
            for (int c = 0; c < variableNodes; ++c) {
                a[r][c].resize(protograph[r][c]);
                for (int id = 0; id < protograph[r][c]; ++id) {
                    a[r][c][id] = getRand(circulant);
                }
            }
        }
        LocalOpt opt(circulant, 6, a);
        int cycles = opt.annealOptimizeNumberOfCycles();


    }
}

int gen(int checkNodes, int variableNodes, const vector<vector<int> >& protograph, ll circulant, ll targetGirth, vector<vector<vector<int> > >& mtr) {
    mtr.assign(checkNodes, vector<vector<int> >(variableNodes));
    for (int r = 0; r < checkNodes; ++r) {
        for (int c = 0; c < variableNodes; ++c) {
            mtr[r][c].resize(protograph[r][c]);
            for (int id = 0; id < protograph[r][c]; ++id) {
                mtr[r][c][id] = getRand(circulant);
            }
        }
    }
    LocalOpt opt(circulant, targetGirth, mtr);
    int cycles = opt.annealOptimizeNumberOfCycles();
    //if (ace != -1) {
    mtr = opt.getMatrix();
    //}
    return cycles;
}

int main(int argc, char* argv[]) {
    /*test();
    return 0;*/

    bool validInput = 1;
    bool regular = 0;
    /*if (argc != 11) {
    validInput = 0;
    }*/
    ll SEED = -1;
    ll GIRTH = -1;
    vector<vector<int> > PROTOGRAPH;
    ll CIRCULANT_SIZE = -1;
    ll VARIABLE_NODES;
    ll CHECK_NODES;
    //ll DESIRED_NUMBER_OF_MATRICES = -1;
    ll NUMBER_OF_ATTEMPTS = -1;

    string INPUT_FILENAME = "";
    for (int i = 1; i + 1 < argc; ++i) {
        if (string(argv[i]) == "-seed") {
            validInput = validInput && toUnsignedInt(argv[i + 1], SEED);
            ++i;
            continue;
        }
        if (string(argv[i]) == "-girth") {
            validInput = validInput && toUnsignedInt(argv[i + 1], GIRTH);
            ++i;
            continue;
        }
        if (string(argv[i]) == "-circulant") {
            validInput = validInput && toUnsignedInt(argv[i + 1], CIRCULANT_SIZE);
            ++i;
            continue;
        }
        if (string(argv[i]) == "-numberOfAttempts") {
            validInput = validInput && toUnsignedInt(argv[i + 1], NUMBER_OF_ATTEMPTS);
            ++i;
            continue;
        }
        if (string(argv[i]) == "-file") {
            INPUT_FILENAME = argv[i + 1];
            ++i;
            continue;
        }
        if (string(argv[i]) == "-regular") {
            validInput = validInput && toUnsignedInt(argv[i + 1], VARIABLE_NODES) && (i + 2 < argc) && toUnsignedInt(argv[i + 2], CHECK_NODES);
            i += 2;
            regular = 1;
            continue;
        }


    }
    if ((GIRTH < 0) || (SEED < 0) || (CIRCULANT_SIZE < 0) || ((INPUT_FILENAME == "") && (!regular)))
        validInput = 0;
    if (!validInput) {
        std::cerr << "Usage: " << argv[0] << " -seed SEED -girth GIRTH -circulant CIRCULANT_SIZE -numberOfMatrices DESIRED_NUMBER_OF_MATRICES -file INPUT_FILENAME" << std::endl;
        return 1;
    }
    if (GIRTH & 1) {
        cerr << "girth must be even\n";
        return 1;
    }

    srand(SEED);
    if (!regular) {
        freopen(INPUT_FILENAME.c_str(), "r", stdin);
        cin >> VARIABLE_NODES >> CHECK_NODES;
        PROTOGRAPH.assign(CHECK_NODES, vector<int>(VARIABLE_NODES));
        for (int i = 0; i < CHECK_NODES; ++i) {
            for (int j = 0; j < VARIABLE_NODES; ++j) {
                cin >> PROTOGRAPH[i][j];
            }
        }
        fclose(stdin);
    }
    else {
        PROTOGRAPH.assign(CHECK_NODES, vector<int>(VARIABLE_NODES, 1));
    }
    if (VARIABLE_NODES % CHECK_NODES != 0) {
        cerr << "WARNING: nonconstant shift\n";

    }
    ll SHIFT = VARIABLE_NODES / CHECK_NODES;
    string folderName = "SpatiallyCoupled" + toStr(VARIABLE_NODES) + "_" + toStr(CHECK_NODES) + "_" + toStr(CIRCULANT_SIZE) + "girth" + toStr(GIRTH);
    string outputFilename = folderName + "/" + toStr(VARIABLE_NODES) + "_" + toStr(CHECK_NODES) + "_" + toStr(CIRCULANT_SIZE) + "girth" + toStr(GIRTH) + "seed" + toStr(SEED);
    if (regular)
        outputFilename += "regular_protograph_matrix";
    else
        outputFilename += "protograph_from_" + INPUT_FILENAME + "_matrix";
    //if (!isPossibleTB(GIRTH, PROTOGRAPH, CIRCULANT_SIZE))
    //return 0;
    system(("mkdir " + folderName).c_str());
    time_t start = time(NULL);
    ll iterationCount = 0;
    //ll successCount = 0;
    ll power = 1;
    int res = -1;
    while ((iterationCount < NUMBER_OF_ATTEMPTS) || (NUMBER_OF_ATTEMPTS == -1)) {
        vector<vector<vector<int> > > a;
        ll cycles = gen(CHECK_NODES, VARIABLE_NODES, PROTOGRAPH, CIRCULANT_SIZE, GIRTH, a);
        if (cycles == -1)
            continue;
        //cerr << "cycles = " << cycles << endl;
        ++iterationCount;
        if ((res == -1) || (res > cycles)) {
            res = cycles;
            freopen((outputFilename + ".txt").c_str(), "w", stdout);
            //++successCount;
            cout << VARIABLE_NODES << "\t" << CHECK_NODES << "\t" << CIRCULANT_SIZE << endl;
            print(a);
            cout << endl;
            fclose(stdout);
            //TB
            //cerr << "girth = " << getGirth(a, PROTOGRAPH, CIRCULANT_SIZE) << endl;
            cerr << "cycles = " << cycles << "\t best = " << res << endl;
            eprint(a);
            cerr << iterationCount << " iterations\n";
            cerr << time(NULL) - start << " seconds\n";
            cerr << endl;
        }
        if (iterationCount == power) {
            power *= 2;
            cerr << time(NULL) - start << " seconds\n";
            cerr << iterationCount << " iterations\n";
            cerr << endl;
        }
    }
    return 0;
}