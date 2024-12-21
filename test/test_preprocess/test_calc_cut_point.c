#include <stdio.h>

#define MAX(A, B) (((A) > (B)) ? (A) : (B))

void find_cut_points(int a[], int n) {
    int cut_point[n];
    int top = 0;
    int remaining_steps = 0;
    for (int i = n - 1; i >= 0; --i) {
        if (remaining_steps == 0) {
            cut_point[top++] = i;
        }
        remaining_steps = MAX(remaining_steps, a[i]+1);
        remaining_steps--;
    }
    int c0 = 0;
    for (int i = top - 1; i >= 0; --i) {
        //printf("%d ", cut_point[i]);
        int c1 = cut_point[i];
        for (int j = c0; j <= c1; ++j) {
            printf("%d ", a[j]);
        }
        printf("\n");
        c0 = c1;
    }
}

int main() {
    int a[] = {0, 0, 2, 0, 1, 1, 3, 2};
    int n = sizeof(a) / sizeof(a[0]);
    find_cut_points(a, n);
    return 0;
}
