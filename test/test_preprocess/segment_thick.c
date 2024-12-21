//
// Created by mainf on 2024/12/1.
//
#include <stdio.h>
#include <stdlib.h>

typedef struct {
    int position;
    int type; // +1 for start, -1 for end
} Event;

int compareEvents(const void *a, const void *b) {
    Event *event1 = (Event *)a;
    Event *event2 = (Event *)b;

    // Sort by position; if positions are the same, end event (-1) comes before start event (+1)
    if (event1->position == event2->position) {
        return event1->type - event2->type;
    }
    return event1->position - event2->position;
}

void calculateThickness(int segments[][2], int num_segments) {
    Event *events = (Event *)malloc(num_segments * 2 * sizeof(Event));

    // Create events for start and end of each segment
    for (int i = 0; i < num_segments; i++) {
        events[i * 2].position = segments[i][0]; // Start
        events[i * 2].type = 1; // +1 for start
        events[i * 2 + 1].position = segments[i][1]; // End
        events[i * 2 + 1].type = -1; // -1 for end
    }

    // Sort events
    qsort(events, num_segments * 2, sizeof(Event), compareEvents);

    int currentThickness = 0;
    int lastPosition = events[0].position;

    printf("Position\tThickness\n");

    // Process events
    for (int i = 0; i < num_segments * 2; i++) {
        // If the position changes, print the thickness
        if (events[i].position != lastPosition) {
            printf("%d\t\t%d\n", lastPosition, currentThickness);
            lastPosition = events[i].position;
        }
        currentThickness += events[i].type; // Update thickness
    }

    // Print the final position thickness
    printf("%d\t\t%d\n", lastPosition, currentThickness);

    // Clean up
    free(events);
}

int main() {
    int segments[][2] = {
        {1, 4},
        {2, 5},
        {3, 6}
    };
    int num_segments = sizeof(segments) / sizeof(segments[0]);

    calculateThickness(segments, num_segments);

    return 0;
}