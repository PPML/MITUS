/////////////////////////// TB REACT PROGRESSION  ////////////////////////////
  for(int ag=0; ag<11; ag++) {
    for(int tb=0; tb<6; tb++) {
      for(int lt=0; lt<2; lt++) {
        for(int nm=0; nm<4; nm++) {
          for(int rg=0; rg<2; rg++) {
            for(int na=0; na<3; na++) {
              temp  = V0[ag][tb][lt][1][nm][rg][na]*vImxtoImyN[ag][0];
              temp2 = V0[ag][tb][lt][2][nm][rg][na]*vImxtoImyN[ag][1];
              temp3 = V0[ag][tb][lt][3][nm][rg][na]*vImxtoImyN[ag][2];
              V1[ag][tb][lt][1][nm][rg][na]  -= temp;
              V1[ag][tb][lt][2][nm][rg][na] += temp;
              V1[ag][tb][lt][2][nm][rg][na] -= temp2;
              V1[ag][tb][lt][3][nm][rg][na]  += temp2;
              V1[ag][tb][lt][3][nm][rg][na] -= temp3;
              V1[ag][tb][lt][4][nm][rg][na]  += temp3;
            } } } } } }
/////////////////////////// NON TB MORTALITY PROGRESSION  ////////////////////
  for(int ag=0; ag<11; ag++) {
    for(int tb=0; tb<6; tb++) {
      for(int lt=0; lt<2; lt++) {
        for(int im=0; im<4; im++) {
          for(int rg=0; rg<2; rg++) {
            for(int na=0; na<3; na++) {
              temp  = V0[ag][tb][lt][im][1][rg][na]*vNmxtoNmyN[ag][0];
              temp2 = V0[ag][tb][lt][im][2][rg][na]*vNmxtoNmyN[ag][1];
              temp3 = V0[ag][tb][lt][im][3][rg][na]*vNmxtoNmyN[ag][2];
              V1[ag][tb][lt][im][1][rg][na]  -= temp;
              V1[ag][tb][lt][im][2][rg][na] += temp;
              V1[ag][tb][lt][im][2][rg][na] -= temp2;
              V1[ag][tb][lt][im][3][rg][na]  += temp2;
              V1[ag][tb][lt][im][3][rg][na] -= temp3;
              V1[ag][tb][lt][im][4][rg][na]  += temp3;
            } } } } } }
