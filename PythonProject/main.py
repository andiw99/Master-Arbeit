# This is a sample Python script.

# Press Shift+F10 to execute it or replace it with your code.
# Press Double Shift to search everywhere for classes, files, tool windows, actions, and settings.
from timeit import default_timer as timer

def fibonacci(n):
    if n < 3:
        return 1
    else:
        return fibonacci(n-1) + fibonacci(n-2)


def main():
    n = int(input("Welche Zahl der Fibonacci Folge möchtest du wissen? "))
    time_before_calc = timer()
    print(f"Die {n}. Fibonacci Zahl lautet {fibonacci(n)}.")
    time_after_calc = timer()
    print(f"Die Berechnung dauerte "
          f"{(time_after_calc - time_before_calc) * 1000:.0f}ms.")

# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    main()

# See PyCharm help at https://www.jetbrains.com/help/pycharm/
