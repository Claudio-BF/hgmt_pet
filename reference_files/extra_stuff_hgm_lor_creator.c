#include <stdio.h>
#include <stdlib.h>
#include "llist.h"

llist* read_eff_line_old(FILE* source) {

    llist* ptr;

    char line[500];
    char* pline;
    double eff;
    char* pEnd;

    // handles if fgets encounters an error
    if (!fgets(line, sizeof(line), source)) {
        return NULL;
    }

    // this points to the beginning character of the line
    pline = line;
    
    while(1){
        
        if (*pline == '\n') {
            break;
        }

        eff = strtod(pline, &pEnd);
        printf("%f\n", eff);
        pline = pEnd;
    }
    
    return ptr;
}

llist* read_eff_line(FILE* source) {

    // variable declaration
    llist* pEffRow = NULL;
    const int entries_in_line = 31;
    int i = 0;
    float eff;
    int worked;

    // loops through all the entries in a row
    while (i < entries_in_line) {

        // assigns the float efficiency found via fscanf to eff
        worked = fscanf(source, "%f,", &eff);

        float* storage = (float*)malloc(sizeof(float));
        *storage = eff;

        // using pointer to the llist, pEffRow, and the pointer to the efficiency, adds element to end of the llist
        pEffRow = add_to_bottom(pEffRow, storage);

        // increments
        i++;
    }
    
    return pEffRow;
}

llist* load_efficiencies(FILE* source, llist* (*f)(FILE*)) {

    llist* pEffTable = NULL;
    llist* pSingleRow = NULL;
    float* eff;
    void* pEff;

    // getting a single row of the table as a llist
    pSingleRow = f(source);
    // making sure my llist is at the top of the whole list
    pSingleRow = list_head(pSingleRow);

    // gets the efficiency of the first element of the llist and prints
    pEff = pSingleRow->data;
    eff = (float*)pEff;
    printf("%f\n", *eff); // this gives the last value in the first row of the .csv
    // printf("%p\n", pSingleRow->data); 

    // gets the efficiency for the second element of the llist (going down) and prints
    llist* next_one_down = pSingleRow->down;
    float* eff_next_down = (float*)(next_one_down->data);
    printf("%f\n", *eff_next_down); // this gives a value that isn't in the .csv at all

    // this segfaults and i've commented it out
    // get_value(pSingleRow, 1);
    eff = (float*)get_value(pSingleRow, 1);
    printf("%f\n", *eff);

    return pEffTable;
}