

#include <iostream>
using namespace std;
struct pol_roots
{
    int num;
    long double* roots;
    int* exp;
};

struct pol
{
    int size;
    long double* data;
};

class Math_solving
{
public:

    struct pol get_pol(long double gain, long double* roots, long double size)
    {
        struct pol res;
        res.data = new long double[2];
        res.size = 1;
        res.data[0] = -gain * roots[0];
        res.data[1] = gain;

        struct pol counter;
        counter.size = 1;
        counter.data = new long double[2];
        for (int i = 1; i < size; i += 1)
        {
            counter.data[0] = -roots[i];
            counter.data[1] = 1;
            struct pol buff = pol_mul(res, counter);

            delete[] res.data;
            res = buff;
        }

        delete[] counter.data;
        return res;
    }

    struct pol pol_sum(struct pol first, struct pol second)
    {
        struct pol res;
        if (first.size > second.size)
            res.size = first.size;
        else
            res.size = second.size;

        res.data = new long double[res.size + 1];
        for (int i = 0; i <= res.size; i += 1)
            res.data[i] = 0;

        for (int i = 0; i <= first.size; i += 1)
            res.data[i] += first.data[i];

        for (int i = 0; i <= second.size; i += 1)
            res.data[i] += second.data[i];

        return res;
    }

    void show_pol(struct pol arg)
    {
        cout << arg.data[0] << " ";
        for (int i = 1; i <= arg.size; i += 1)
            cout << (arg.data[i] >= 0 ? "+ " : "- ") << (arg.data[i] < 0 ? arg.data[i] * -1 :
                arg.data[i]) << " * X^" << i << " ";
    }

    struct pol pol_mul(struct pol first, struct pol second)
    {
        struct pol res;
        res.size = first.size + second.size;
        res.data = new long double[res.size + 1];
        for (int i = 0; i <= res.size; i += 1)
            res.data[i] = 0;

        for (int i = 0; i <= second.size; i += 1)
            for (int i2 = 0; i2 <= first.size; i2 += 1)
                res.data[i + i2] += second.data[i] * first.data[i2];

        return res;
    }

    struct pol_roots real_solve_pol(struct pol arg, long double toler)
    {
        int i = arg.size;
        while (i > 0 and arg.data[i] == 0)
            i -= 1;

        if (i < 0)
            cout << "Empty pol!\n";
        else
        {
            arg.size = i;
            for (int i = 0; i <= arg.size; i += 1)
                arg.data[i] /= arg.data[arg.size];

            return real_solve_pol_procedure(arg.data, arg.size, toler);
        }
    }

    long double get_value_pol(long double X, long double* pol, int size)
    {
        long double Gain = 1;
        long double Res = 0;

        for (int i = 0; i <= size; i += 1)
        {
            Res += Gain * pol[i];
            Gain *= X;
        }

        return Res;
    }

private:

    struct pol_roots real_solve_pol_procedure(long double* pol, int size, long double toler)
    {
        if (size > 2)
        {
            int dif_size = size - 1;
            long double* dif_pol = new long double[dif_size + 1];

            for (int i = 1; i <= size; i += 1)
                dif_pol[i - 1] = pol[i] * i;

            struct pol_roots dif_roots = real_solve_pol_procedure(dif_pol, dif_size, toler);
            delete[] dif_pol;

            int extr_points_num = dif_roots.num;
            long double* extr_points_val = new long double[extr_points_num + 2];
            long double* extr_points = new long double[extr_points_num + 2];
            int* extr_points_exp = new int[extr_points_num + 2];
            int start;

            long double sign;
            if ((pol[size] > 0 and size % 2 == 0) or (pol[size] < 0 and size % 2 == 1))
                sign = 1;
            else
                sign = -1;

            long double left;
            if (dif_roots.num > 0)
                left = get_value_pol(dif_roots.roots[0], pol, size);

            if ((dif_roots.num > 0 and sign * left < 0) or dif_roots.num == 0)
            {
                long double buff, buff_value;
                long double area = 1;
                if (dif_roots.num == 0)
                    buff = 0;
                else
                    buff = dif_roots.roots[0] - area;

                buff_value = get_value_pol(buff, pol, size);
                while (left * buff_value > 0)
                {
                    area *= 4;
                    buff -= area;
                    buff_value = get_value_pol(buff, pol, size);
                }

                extr_points[0] = buff;
                extr_points_val[0] = buff_value;
                extr_points_exp[0] = 0;
                extr_points_num += 1;
                start = 1;
            }
            else
                start = 0;

            for (int i = 0; i < dif_roots.num; i += 1)
            {
                extr_points_val[i + start] = get_value_pol(dif_roots.roots[i], pol, size);
                extr_points[i + start] = dif_roots.roots[i];
                extr_points_exp[i + start] = dif_roots.exp[i];
            }

            if ((pol[size] > 0 and size % 2 == 0) or (pol[size] > 0 and size % 2 == 1))
                sign = 1;
            else
                sign = -1;

            if (sign * extr_points_val[start + dif_roots.num - 1] < 0)
            {
                long double buff, buff_value;
                long double area = 1;
                buff = extr_points[start + dif_roots.num - 1] + area;

                buff_value = get_value_pol(buff, pol, size);
                while (extr_points_val[start + dif_roots.num - 1] * buff_value > 0)
                {
                    area *= 4;
                    buff += area;
                    buff_value = get_value_pol(buff, pol, size);
                }

                extr_points[start + dif_roots.num] = buff;
                extr_points_val[start + dif_roots.num] = buff_value;
                extr_points_exp[start + dif_roots.num] = 0;
                extr_points_num += 1;
            }

            struct pol_roots res;
            res.num = 0;
            res.roots = new long double[size];
            res.exp = new int[size];

            for (int i = 0; i < extr_points_num; i += 1)
            {
                if (-toler <= extr_points_val[i] and extr_points_val[i] <= toler)
                {
                    res.roots[res.num] = extr_points[i];
                    res.exp[res.num] = extr_points_exp[i] + 1;
                    res.num += 1;
                }
                else  if ((-toler > extr_points_val[i + 1] or extr_points_val[i + 1] > toler) and i < (extr_points_num - 1))
                {
                    if ((extr_points_val[i] < 0 and extr_points_val[i + 1] > 0) or
                        (extr_points_val[i] > 0 and extr_points_val[i + 1] < 0))
                    {
                        res.roots[res.num] = binary_iter_pol(extr_points[i], extr_points[i + 1], pol, size, toler);
                        res.exp[res.num] = 1;
                        res.num += 1;
                    }
                }
            }

            delete[] dif_roots.roots;
            delete[] dif_roots.exp;
            delete[] extr_points_exp;
            delete[] extr_points_val;
            delete[] extr_points;
            return  res;
        }
        else if (size == 2)
        {
            struct pol_roots res;
            res.roots = new long double[2];
            res.exp = new int[2];
            long double D = pol[1] * pol[1] - 4 * pol[2] * pol[0];

            if (D < 0)
                res.num = 0;
            else
            {
                long double sqrt_D = sqrt(D);
                res.roots[0] = (-pol[1] - sqrt_D) / (2 * pol[2]);
                res.num = 1;
                if (D != 0)
                {
                    res.roots[1] = (-pol[1] + sqrt_D) / (2 * pol[2]);
                    res.num = 2;
                    res.exp[0] = 1;
                    res.exp[1] = 1;
                }
                else
                    res.exp[0] = 2;
            }

            return res;
        }
        else
        {
            struct pol_roots res;
            res.num = 1;
            res.roots = new long double[1];
            res.exp = new int[1];

            res.roots[0] = -pol[0] / pol[1];
            res.exp[0] = 1;
            return  res;
        }
    }

    long double binary_iter_pol(long double A, long double B, long double* pol, int size, long double toler)
    {
        long double C, pol_A, pol_B, pol_C;

        while (1)
        {
            pol_A = get_value_pol(A, pol, size);
            pol_B = get_value_pol(B, pol, size);

            C = (B + A) / 2;
            pol_C = get_value_pol(C, pol, size);

            if (-toler <= pol_C and pol_C <= toler)
                return C;
            else
            {
                if (pol_A * pol_C < 0)
                    B = C;
                else
                    A = C;
            }
        }
    }
};

int main()
{
    class Math_solving test;
    long double toler = pow(10, -7);

    int size = 6;
    long double* roots = new long double[size];
    long double gain = 85;

    roots[0] = -14;
    roots[1] = -14;
    roots[2] = -14.7;
    roots[3] = -14.7;
    roots[4] = -14;
    roots[5] = -14;

    struct pol arg = test.get_pol(gain, roots, size);
    delete[] roots;
    test.show_pol(arg);

    struct pol_roots res = test.real_solve_pol(arg, toler);
    cout << "\n\nroots:\n";
    for (int i = 0; i < res.num; i += 1)
        cout << res.roots[i] << "    exp: " << res.exp[i] << endl;
    cout << endl;
    delete[] res.roots;

}

