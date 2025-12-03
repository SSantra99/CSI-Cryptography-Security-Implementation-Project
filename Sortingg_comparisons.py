import matplotlib.pyplot as plt
import numpy as np

n = np.arange(1, 2001, 50)  # Array sizes

# Hypothetical time complexities (normalized)
bubble = n**2
selection = n**2
insertion = n**2
merge = n * np.log2(n)
quick = n * np.log2(n)
heap = n * np.log2(n)
counting = n  # assuming k ~ n
radix = n    # assuming k ~ n

plt.figure(figsize=(10,6))
plt.plot(n, bubble, label='Bubble Sort')
plt.plot(n, selection, label='Selection Sort')
plt.plot(n, insertion, label='Insertion Sort')
plt.plot(n, merge, label='Merge Sort')
plt.plot(n, quick, label='Quick Sort')
plt.plot(n, heap, label='Heap Sort')
plt.plot(n, counting, label='Counting Sort')
plt.plot(n, radix, label='Radix Sort')

plt.xlabel('Array Size (n)')
plt.ylabel('Time Complexity (normalized)')
plt.title('Sorting Algorithms Time Complexity Comparison')
plt.legend()
plt.grid(True)
plt.show()
