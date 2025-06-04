#include <stdio.h>
#include <stdlib.h>
#include "llist.h"

int main() {

    // define an integer i and sets value to 5
    int i;
    i = 5;
    // prints result
    printf("Value of i: %i\n", i);

    // define a pointer pi and sets to address of i
    int *pi;
    pi = &i;
    printf("Value of i (using pointer): %i\n", *pi);
    printf("Address of i: %p\n", pi);

    // changes the value of i using the pointer
    *pi = 6;
    printf("Value of i (changed using pointer): %i\n", i);
    printf("Address of i: %p\n", pi); 

    // creates a pointer called j
    int *j = (int*)malloc(sizeof(int));
    *j = 4;
    printf("Value of pointer j: %p\n", j);
    int* k = (int*)calloc(j[0], sizeof(int));
    for (int i = 0; i < j[0]; i++) {
        k[i] = i*i;
    }
    for (int i = j[0] - 1; i >= 0; i--) {
        printf("%i", k[i]);
    }
    free(j);




    // creates pointer called my_llist to a llist object
    // my_llist is the address of the object in memory
    // *my_llist is the object itself 
    // llist *my_llist;

    //add_to_bottom(my_llist, 5);

    return 0;
}