/**
 *
 * AUTHOR - Parthasarathi Das
 * BRIEF  - A class for quadratic ideals
 */

#ifndef QF_HPP
#define QF_HPP

#include <NTL/ZZ_limbs.h>

NTL_CLIENT

class QF
{
    public:
    QF();
    QF(mpz_qform_group_t& group);
    ~QF();
    
    void init            (mpz_qform_group_t& group);
    void assign          (const mpz_t a, const mpz_t b, const mpz_t D);
    void getb            (mpz_t b);
    void setD            (const mpz_t D);
    void getD            (mpz_t D);
    void random_prime_QF ();
    int  next_prime_QF   (int& index, const mpz_t D);
    void findQFprimeto   (const QF& I, const mpz_t D, const mpz_t f);
    void clear           (mpz_qform_group_t& group);
    
    friend void assign   (QF& Z, const QF& X);
    friend void clear    (QF& Z);
    friend bool equals   (QF& Z, const QF& X);
    friend void inv      (QF& Z, const QF& X);
    friend void mul      (QF& Z, const QF& X, const QF& Y); // forms composition
    friend void sqr      (QF& Z, const QF& X);
    friend void cub      (QF& Z, const QF& X);
    friend void liftQF   (QF& Z, const QF& I, const mpz_t DK, const mpz_t D, const mpz_t f);
    
    friend std::ostream& operator<< (std::ostream& out, const QF& A);
    friend bool          operator== (const QF& X, const QF& Y);
    
    
    protected:
    mpz_qform_group_t* gp;
    mpz_qform_t form;
};

// Out of line class member definitions
inline QF::QF ()
{
    mpz_qform_init (gp, &form);
}

inline QF::~QF ()
{
    mpz_qform_clear (gp, &form);
}

inline QF::QF (mpz_qform_group_t& group)
{
    gp = &group;
    mpz_qform_init (gp, &form);
}

inline void QF::init (mpz_qform_group_t& group)
{
    gp = &group;
}

inline void QF::clear (mpz_qform_group_t& group)
{
    mpz_qform_group_clear (&group);
}

inline void QF::assign (const mpz_t a, const mpz_t b, const mpz_t D)
{
    mpz_qform_group_set_discriminant (gp, D);
    mpz_set (form.a, a);
    mpz_set (form.b, b);
    mpz_qform_c (gp, form.c, form.a, form.b);
}

inline void QF::getb (mpz_t b)
{
    mpz_set (b, form.b);
}
inline void QF::setD (const mpz_t D)
{
    mpz_qform_group_set_discriminant (gp, D);
}

inline void QF::getD (mpz_t D)
{
    mpz_set (D, gp->D);
}

inline void QF::random_prime_QF ()
{
    qform_random_primeform (&gp->desc, &form);
}

inline int QF::next_prime_QF (int& index, const mpz_t D)
{
    mpz_qform_group_set_discriminant (gp, D);
    return qform_next_primeform (&gp->desc, &form, index);
}

inline void QF::findQFprimeto(const QF& I, const mpz_t D, const mpz_t f)
{    
    mpz_qform_group_set_discriminant (gp, D);
    mpz_gcd (form.c, I.form.a, f); // using form.c as a buffer
    
    if(mpz_cmp_ui (form.c, 1) > 0)
    {
        mpz_gcd (form.c, I.form.c, f);
        if(mpz_cmp_ui(form.c, 1) > 0)
        {
            mpz_add (form.a, I.form.a, I.form.b);
            mpz_add (form.a, form.a, I.form.c);
            mpz_mul_ui (form.c, I.form.a, 2);
            mpz_add (form.b, I.form.b, form.c);
            mpz_neg (form.b, form.b);
        }
        else
        {
            mpz_set (form.a, I.form.c);
            mpz_set (form.b, I.form.b);
            mpz_neg (form.b, form.b);
        }
        
        mpz_qform_c (gp, form.c, form.a, form.b); // here, form.c is calculated
    }
    else
    {
        mpz_set (form.a,I.form.a);
        mpz_set (form.b,I.form.b);
        mpz_qform_c (gp, form.c, form.a, form.b); // here, form.c is calculated
    }
}

// Non class members
inline void assign (QF& Z, const QF& X)
{
    Z.gp = X.gp;
    mpz_qform_set (Z.gp, &Z.form, &X.form);
}

inline void clear (QF& Z)
{
    mpz_qform_clear (Z.gp, &Z.form);
}

inline bool equals (QF& Z, const QF& X)
{
    return mpz_qform_equal (Z.gp, &Z.form, &X.form);
}

inline void inv (QF& Z, const QF& X)
{
    assign(Z, X);
    mpz_qform_inverse (Z.gp, &Z.form);
}

inline void mul (QF& Z, const QF& X, const QF& Y)
{
    Z.gp = X.gp;
    mpz_qform_compose (Z.gp, &Z.form, &X.form, &Y.form);
}

inline void sqr (QF& Z, const QF& X)
{
    Z.gp = X.gp;
    mpz_qform_square (Z.gp, &Z.form, &X.form);
}

inline void cub (QF& Z, const QF& X)
{
    Z.gp = X.gp;
    mpz_qform_cube (Z.gp, &Z.form, &X.form);
}

inline void liftQF(QF& Z, const QF& I, const mpz_t Dsmall, const mpz_t Dlarge, const mpz_t f)
{
    Z.findQFprimeto (I, Dsmall, f); // Here group has D = Dsmall;
    mpz_mul (Z.form.b, Z.form.b, f);
    mpz_mul_2exp (Z.form.c, Z.form.a, 1); // Using Z.form.c as a buffer
    mpz_mod (Z.form.b, Z.form.b, Z.form.c); // Using Z.form.c as a buffer
    mpz_qform_group_set_discriminant (Z.gp, Dlarge); // Now group has D = Dlarge
    mpz_qform_c (Z.gp, Z.form.c, Z.form.a, Z.form.b); // Z.form.c is calculated
}

inline std::ostream& operator<< (std::ostream& out, const QF& A)
{
    mpz_qform_print (A.gp, &A.form);
    gmp_printf("\n");
    return out;
}

inline bool operator== (const QF& X, const QF& Y)
{
    return mpz_qform_equal (X.gp, &X.form, &Y.form);
}

#endif // QF_HPP


