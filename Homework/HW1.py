import numpy as np
import matplotlib.pyplot as plt

y1 = lambda x,d: np.cos(x+d) - np.cos(x)

y2 = lambda x, d: -2*np.sin(x + d/2)*np.sin(d/2)

delta = []
for i in range(-16,1):
    delta.append(10**i)

# x = 2pi
x = np.pi*2
err = []
for d in delta:
    err.append(abs(y1(x,d)-y2(x,d)))

print(delta)
print(err)

plt.figure(1)
plt.plot(delta, err)
plt.xscale('log')
plt.title('$x=2\pi$')
plt.ylabel('error')
plt.ylabel('\delta')
plt.show()

# x = 10^6
x = 10**6
err = []
for d in delta:
    err.append(abs(y1(x,d)-y2(x,d)))

print(delta)
print(err)

plt.figure(2)
plt.plot(delta, err)
plt.xscale('log')
plt.title('$x=10^{6}$')
plt.ylabel('error')
plt.ylabel('\delta')
plt.show()
