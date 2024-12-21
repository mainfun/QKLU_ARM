#include <stdio.h>
#include <stdlib.h>

typedef struct {
    int value;
    int arrayIndex; // 表示该元素来自哪个有序段
    int elementIndex; // 表示该元素在其段内的索引
} HeapNode;

// 最小堆的实现
typedef struct {
    HeapNode *nodes;
    int size;
} MinHeap;

MinHeap* createMinHeap(int capacity) {
    MinHeap *heap = (MinHeap *)malloc(sizeof(MinHeap));
    heap->nodes = (HeapNode *)malloc(capacity * sizeof(HeapNode));
    heap->size = 0;
    return heap;
}

void swap(HeapNode *a, HeapNode *b) {
    HeapNode temp = *a;
    *a = *b;
    *b = temp;
}

void minHeapify(MinHeap *heap, int index) {
    int smallest = index;
    int left = 2 * index + 1;
    int right = 2 * index + 2;

    if (left < heap->size && heap->nodes[left].value < heap->nodes[smallest].value) {
        smallest = left;
    }
    if (right < heap->size && heap->nodes[right].value < heap->nodes[smallest].value) {
        smallest = right;
    }
    if (smallest != index) {
        swap(&heap->nodes[index], &heap->nodes[smallest]);
        minHeapify(heap, smallest);
    }
}

void insertMinHeap(MinHeap *heap, HeapNode node) {
    heap->nodes[heap->size] = node;
    heap->size++;
    for (int i = (heap->size - 1) / 2; i >= 0; i--) {
        minHeapify(heap, i);
    }
}

HeapNode extractMin(MinHeap *heap) {
    HeapNode root = heap->nodes[0];
    heap->nodes[0] = heap->nodes[--heap->size];
    minHeapify(heap, 0);
    return root;
}

int isEmpty(MinHeap *heap) {
    return heap->size == 0;
}

// 归并多个有序段
void mergeKSortedArrays(int **arrays, int *sizes, int k, int *output, int *outputSize) {
    MinHeap *minHeap = createMinHeap(k);

    // 初始化堆
    for (int i = 0; i < k; i++) {
        if (sizes[i] > 0) {
            HeapNode node;
            node.value = arrays[i][0];
            node.arrayIndex = i;
            node.elementIndex = 0;
            insertMinHeap(minHeap, node);
        }
    }

    *outputSize = 0;

    while (!isEmpty(minHeap)) {
        HeapNode minNode = extractMin(minHeap);
        output[(*outputSize)++] = minNode.value;

        // 如果该节点的数组中还有下一个元素，则将其插入堆中
        if (minNode.elementIndex + 1 < sizes[minNode.arrayIndex]) {
            HeapNode nextNode;
            nextNode.value = arrays[minNode.arrayIndex][minNode.elementIndex + 1];
            nextNode.arrayIndex = minNode.arrayIndex;
            nextNode.elementIndex = minNode.elementIndex + 1;
            insertMinHeap(minHeap, nextNode);
        }
    }

    free(minHeap->nodes);
    free(minHeap);
}

int main() {
    int k = 3; // 有序段的数量
    int sizes[] = {3, 3, 3}; // 每个段的大小
    int *arrays[3];

    // 初始化有序段
    int arr1[] = {1, 4, 7};
    int arr2[] = {2, 5, 8};
    int arr3[] = {3, 6, 9};
    arrays[0] = arr1;
    arrays[1] = arr2;
    arrays[2] = arr3;

    int output[9]; // 输出数组
    int outputSize = 0;

    mergeKSortedArrays(arrays, sizes, k, output, &outputSize);

    // 打印结果
    printf("Merged sorted array:\n");
    for (int i = 0; i < outputSize; i++) {
        printf("%d ", output[i]);
    }
    printf("\n");

    return 0;
}