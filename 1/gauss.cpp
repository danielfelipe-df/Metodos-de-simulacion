#include <iostream>
#include <fstream>
#include <string>
#include <vector>

int main(void)
{
    char *locale = setlocale(LC_ALL, "");
    char c;

    FILE * pFile2 = fopen("datos_5g.dat", "r");
    std::vector<std::string> w;
    int k=0;
    while ((c = fgetc(pFile2)) != WEOF){
      if(iswalnum(c)){
        w[k].push_back(c);
      }
      else{
        k++;
        w[k].push_back(c);
        k++;
      }
    }
    fclose(pFile2);

    for(int i=0; i<w.size(); i++){
        std::cout << w[i] << std::endl;
    }

    return 0;
}
