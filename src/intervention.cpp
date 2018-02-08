//////////////////////// INTERVENTION ENROLLMENT /////////////////////////////
// rate differs by RISK category, rg group
for(int ag=0; ag<11; ag++) {
    for(int tb=0; tb<6; tb++) {
        for(int lt=0; lt<2; lt++) {
            for(int rg=0; rg<2; rg++) {
                for(int na=0; na<3; na++) {
                    for(int rf2=0; rf2<4; rf2++) {
                        if( tb>5 | rf2==1 ) {
                            ti = 1;
                        } else { ti = 0;
                        }   // ti = 0 for slow rate, 1 for higher rate (advanced HIV and/or diag TB)
                        if( rg==1 ) {
                            ti += 2; } // 3rd and 4th col of rIntvInit for HR group (will need to update)
                        temp = V0[ag][tb][lt][1+rf2][1+rf2][rg][na]*rIntvInitN[s][ti];
                        V1[ag][tb][lt][1+rf2][1+rf2][rg][na]  -= temp;
                        V1[ag][tb][lt][2+rf2][2+rf2][rg][na]  += temp;
                    } } } } } }
////////////////////////  INTERVENTION DEFAULT  /////////////////////////////
for(int ag=0; ag<11; ag++) {
    for(int tb=0; tb<6; tb++) {
        for(int lt=0; lt<2; lt++) {
            for(int rg=0; rg<2; rg++) {
                for(int na=0; na<3; na++) {
                    if(rg!=1) { temp3=rIntvDef; } else { temp3=rIntvDef*2; }
                    temp  = V0[ag][tb][dr][tx][2 ][rg]*temp3;
                    temp2 = V0[ag][tb][dr][tx][4 ][rg]*temp3;
                    V1[ag][tb][dr][tx][2 ][rg]  -= temp;
                    V1[ag][tb][dr][tx][1 ][rg]  += temp;
                    V1[ag][tb][dr][tx][4 ][rg]  -= temp2;
                    V1[ag][tb][dr][tx][3 ][rg]  += temp2;
                } } } }

