//
// Created by mainf on 2024/12/18.
//
// 自定义比较器，基于 values[i] 进行比较
struct Compare {
    const INDEX_TYPE *values;

    Compare(const INDEX_TYPE *vals) : values(vals) {
    }

    // 小值优先
    bool operator()(INDEX_TYPE a, INDEX_TYPE b) {
        return values[a] > values[b]; // 如果 values[a] 大于 values[b]，则 a 优先级低
    }
};

INDEX_TYPE *toposort_by_etree_v3(INDEX_TYPE n, const INDEX_TYPE parent[],
                                 const INDEX_TYPE values[], const INDEX_TYPE threshold,
                                 INDEX_TYPE *cut_point, const INDEX_TYPE values2[]) {
    INDEX_TYPE *order = (INDEX_TYPE *) lu_malloc(n * sizeof(INDEX_TYPE));
    std::stack<INDEX_TYPE> small_stack; // 大值优先队列
    std::stack<INDEX_TYPE> big_queue;
    INDEX_TYPE *kid_count = (INDEX_TYPE *) lu_calloc(n + 1, sizeof(INDEX_TYPE));
    INDEX_TYPE top_order = -1;
    kid_count++;

    // 计算每个节点的子节点数
    for (INDEX_TYPE i = 0; i < n; kid_count[parent[i++]]++);

    // 将所有叶子节点入栈或入队
    for (INDEX_TYPE i = n - 1; i >= 0; i--) {
        if (kid_count[i] == 0) { // 叶子节点
            if (values[i] > threshold) {
                big_queue.push(i); // 加入大值队列
            } else {
                small_stack.push(i); // 入小值栈
            }
        }
    }

    bool is_first_pop_big_queue = true;

    // 处理队列中的节点
    while (!small_stack.empty() || !big_queue.empty()) {
        INDEX_TYPE node;

        // 从小值栈取出一个节点，若空则从大值队列取
        if (!small_stack.empty()) {
            node = small_stack.top();
            small_stack.pop(); // 弹出小值栈
        } else {
            node = big_queue.top();
            big_queue.pop(); // 弹出大值队列
            if (is_first_pop_big_queue) {
                *cut_point = node;
                is_first_pop_big_queue = false;
            }
        }

        order[++top_order] = node;

        // 获取父节点，并减少其子节点数
        INDEX_TYPE parent_node = parent[node];
        kid_count[parent_node]--;

        // 如果父节点的子节点数为0，检查其值并入队
        if (parent_node > 0 && kid_count[parent_node] == 0) {
            if (values[parent_node] > threshold) {
                big_queue.push(parent_node); // 加入大值队列
            } else {
                small_stack.push(parent_node); // 入小值栈
            }
        }
    }
    kid_count--;
    delete kid_count;
    return order;
}

//只排小值栈
INDEX_TYPE *toposort_0(INDEX_TYPE n, const INDEX_TYPE parent[]
                       // const INDEX_TYPE values[], const INDEX_TYPE threshold,
                       // INDEX_TYPE *cut_point
) {
    INDEX_TYPE *order = (INDEX_TYPE *) lu_malloc(n * sizeof(INDEX_TYPE));
    INDEX_TYPE *sibling_sort = (INDEX_TYPE *) lu_malloc(n * sizeof(INDEX_TYPE));
    INDEX_TYPE sort_top = 0, order_top = 0;
    INDEX_TYPE *first_kid = (INDEX_TYPE *) lu_malloc((n + 1) * sizeof(INDEX_TYPE));
    INDEX_TYPE *sibling = (INDEX_TYPE *) lu_malloc((n + 1) * sizeof(INDEX_TYPE));
    for (INDEX_TYPE v = 0; v <= n; v++) {
        first_kid[v] = sibling[v] = -1;
    }
    first_kid++, sibling++;
    parent2child_sibling(parent, n, first_kid, sibling);
    INDEX_TYPE current = 0;
    while (current != -1) {
        sort_top = 0;
        INDEX_TYPE node = current;
        do {
            sibling_sort[sort_top++] = node; //todo remove sibling_sort
            node = sibling[node];
        } while (node != -1);
        std::sort(sibling_sort, sibling_sort + sort_top);
        for (INDEX_TYPE i = 0; i < sort_top; i++) {
            INDEX_TYPE v = sibling_sort[i];
            order[order_top++] = v;
        }
        current = parent[current];
    }

    first_kid--, sibling--;
    lu_free(sibling_sort);
    lu_free(sibling);
    lu_free(first_kid);
    return order;
}


///父亲表示法转孩子兄弟表示法
void parent2child_sibling(const INDEX_TYPE *parent, const INDEX_TYPE n,
                          INDEX_TYPE *first_kid,INDEX_TYPE *sibling) {
    for (INDEX_TYPE v = 0; v < n; v++) {
        printf("%lld,", parent[v]);
        INDEX_TYPE dad = parent[v];
        sibling[v] = first_kid[dad];
        first_kid[dad] = v;
    }
}