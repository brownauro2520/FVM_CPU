void residual(vector<double> &R, vector<double> &dT, vector<double> U, vector<vector<double>> I2E, vector<vector<double>> B2E, vector<vector<double>> In, vector<vector<double>> Bn,
	vector<double> Il, vector<double> Bl, vector<double> A, vector<double> u_inf, int N, double CFL) {

    double gamma = 1.4; double gmi = gamma - 1;

    int ul; int ur; double dl;
    vector<double> UL(4);
    vector<double> UR(4);
    vector<double> n(2);
    pair<vector<double>, double> FSI({ 0.0, 0.0, 0.0, 0.0 }, 0.0);

    for (int i = 0; i < I2E.size(); i++) {

        //obtaining the index for I2E
        ul = I2E[i][0] - 1;
        ur = I2E[i][2] - 1;

        //obtaining the states of the elmemnts
        UL = { U[ul * 4], U[ul * 4 + 1], U[ul * 4 + 2], U[ul * 4 + 3] };
        UR = { U[ur * 4], U[ur * 4 + 1], U[ur * 4 + 2], U[ur * 4 + 3] };

        //obtaiining normal
        n = In[i];

        //caling flux function
        FSI = HLLE(UL, UR, n, gamma);

        //obtaiing length
        dl = Il[i];

        //adding the to the element
        for (int j = 0; j < 4; j++) {
            R[ul * 4 + j] += FSI.first[j] * dl;
        }

        //subtaracting from the element
        for (int j = 0; j < 4; j++) {
            R[ur * 4 + j] -= FSI.first[j] * dl;
        }

        //computing dT
        dT[ul] += abs(FSI.second) * dl;
        dT[ur] += abs(FSI.second) * dl;

    }

    //--------------------------------------------------------------------------------------------------------------------------
    //boundary edges resiudal computation

    int check;

    for (int i = 0; i < B2E.size(); i++) {

        //ibtaining elmenet index
        ul = B2E[i][0] - 1;

        //obtaining boundary group
        check = B2E[i][2];

        //obtaining states of the element
        UL = { U[ul * 4], U[ul * 4 + 1], U[ul * 4 + 2], U[ul * 4 + 3] };

        //obtaing normal
        n = Bn[i];

        //computing flux
        //farfield - full state condition
        if (check == -1) {
            FSI = HLLE(UL, u_inf, n, gamma);

        }

        //airfoil boundary
        else {
            FSI = inviscid(UL, n);
        }

        //obtaiing length
        dl = Bl[i];

        // adding the to the element
        for (int j = 0; j < 4; j++) {
            R[ul * 4 + j] += FSI.first[j] * dl;
        }

        //computing dT
        dT[ul] += abs(FSI.second) * dl;
    }

    //--------------------------------------------------------------------------------------------------------------------------
    //computing time stepping
    for (int i = 0; i < N; i++) {
        dT[i] = (2 * CFL * A[i]) / dT[i];
    }
}