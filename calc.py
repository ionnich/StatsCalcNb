from scipy import stats

def tTable(df, alpha, tail):
    if(tail == "two"):
        return stats.t.ppf(1-alpha/2, df)

    elif(tail == "left"):
        return stats.t.ppf(alpha, df)

    elif(tail == "right"):
        return stats.t.ppf(1-alpha, df)

def zTable(alpha, tail):
    if(tail == "two"):
        return stats.norm.ppf(1-alpha/2)

    elif(tail == "left"):
        return stats.norm.ppf(alpha)

    elif(tail == "right"):
        return stats.norm.ppf(1-alpha)

def hypothesisTwoPopulationsZ(x1, x2, s1, s2, n1, n2, alpha, tail):
    # x1 = mean of sample 1
    # x2 = mean of sample 2
    # s1 = std of sample 1
    # s2 = std of sample 2
    # n1 = size of sample 1
    # n2 = size of sample 2
    # alpha = significance level
    # tail = two, left or right
    # return z, p, reject

    # z statistic
    z = (x1 - x2) / math.sqrt(s1**2/n1 + s2**2/n2)

    # critical value
    c = zTable(alpha, tail)

    if(tail == "two"):
        if(z > zTable(alpha/2, "right")):
            reject = True
            c = zTable(alpha/2, "right")
        if(z < zTable(alpha/2, "left")):
            reject = True
            c = zTable(alpha/2, "left")
        else:
            reject = False

    elif(tail == "left"):
        if(z < zTable(alpha, "left")):
            c = zTable(alpha, "left")
            reject = True

    elif(tail == "right"):
        if(z > zTable(alpha, "right")):
            c = zTable(alpha, "right")
            reject = True

    return z, c, reject

# for t-table testing
def hypothesisTwoPopulationsT(x1, x2, s1, s2, n1, n2, alpha, tail):
    # x1 = mean of sample 1
    # x2 = mean of sample 2
    # s1 = std of sample 1
    # s2 = std of sample 2
    # n1 = size of sample 1
    # n2 = size of sample 2
    # alpha = significance level
    # tail = two, left or right
    # return t, p, reject

    # t statistic
    t = (x1-x2) / (((s1**2)/n1) + ((s2**2)/n2))**(1/2)

    # ask if population standard deviations are equal
    equal = input("Are the population standard deviations equal? (y/n): ")

    df = n1 + n2 - 2
    if(equal == "y"):
        # pooled standard deviation
        sp = ((n1-1)*s1**2 + (n2-1)*s2**2) / (n1+n2-2)
        sp = sp**(1/2)
        # print this formula in latex
        print("t = \\frac{\\bar{x}_1 - \\bar{x}_2}{\\sqrt{\\frac{s_1^2}{n_1} + \\frac{s_2^2}{n_2}}}")

        # t statistic
        t = (x1-x2) / (sp*((1/n1)+(1/n2)))**(1/2)
        # degrees of freedom

    # critical value
    c = tTable(df, alpha, tail)
    
    reject = False

    if(tail == "two"):
        if(t > c):
            reject = True
        if(t < -c):
            c = -c 
            reject = True
    if(tail == "left"):
        if(t < -c):
            c = -c
            reject = True
    if(tail == "right"):
        if(t > c):
            reject = True

    return t, c, reject

def hypothesisTwoPopulations(x1, x2, s1, s2, n1, n2, alpha, tail):
    # x1 = mean of sample 1
    # x2 = mean of sample 2
    # s1 = std of sample 1
    # s2 = std of sample 2
    # n1 = size of sample 1
    # n2 = size of sample 2
    # alpha = significance level
    # tail = two, left or right
    # return t, p, reject

    if(n1 > 30 and n2 > 30):
        return hypothesisTwoPopulationsZ(x1, x2, s1, s2, n1, n2, alpha, tail)

    else:
        return hypothesisTwoPopulationsT(x1, x2, s1, s2, n1, n2, alpha, tail)



if __name__ == "__main__":

    print("Welcome to the Hypothesis Testing Calculator!")
    print("This program will calculate the t-statistic and p-value for a hypothesis test.")


    # get sample data
    x1 = float(input("Enter the mean of sample 1: "))
    x2 = float(input("Enter the mean of sample 2: "))
    s1 = float(input("Enter the standard deviation of sample 1: "))
    s2 = float(input("Enter the standard deviation of sample 2: "))
    n1 = int(input("Enter the size of sample 1: "))
    n2 = int(input("Enter the size of sample 2: "))
    alpha = float(input("Enter the significance level (alpha): "))
    tail = input("Enter the tail (two, left, right): ")

    # calculate t-statistic and p-value
    t, c, reject = hypothesisTwoPopulations(x1, x2, s1, s2, n1, n2, alpha, tail)

    # print results
    if(n1 > 30 and n2 > 30):
        print("z = ", t)
    else:
        print("t = ", t)
    print("critical value: ", c)
    print("reject: ", reject)

    print("Would you like to input new data? (y/n)")
    again = input()
    if(again == "y"):
        main()

    pass
