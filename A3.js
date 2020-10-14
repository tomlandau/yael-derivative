const mathjs = require('mathjs');
const Calculess = require('calculess');
const Calc = Calculess.prototype;

let Xb = 1;
let Ya = 1;
let Yb = 0;
let Za = 0;
let Zb = 1;

let h1 = 2.925;
let h2 = 0.325;
let n1 = 2; // Quanti
let n2 = 1; // Quanti
let l1 = 2; // Quanti
let l2 = 0; // Quanti
let lamda = 0; // positive, same as m1
let m2 = 0;

const delta = (lamda != m2) ? 0 : 1;
const t = h2 / h1;
const N = n1 + n2 + 1;

const H = () => {
    const divider = mathjs.sqrt(mathjs.factorial(2*n1) * mathjs.factorial(2 * n2));
    const firstPow = Math.pow(-1, l2-lamda);
    const secondPow =  Math.pow(1+t,n1+0.5);
    const thirdPow = Math.pow(1-t, n2+0.5);

    return (delta * firstPow * secondPow * thirdPow) / divider;
};

const R = (Xa) => {
    return mathjs.sqrt(Math.pow(Xa-Xb, 2) + Math.pow(Ya-Yb, 2) + Math.pow(Za-Zb, 2));
}

const ro = (Xa) => {
    return R(Xa) * h1/2;
};

const F = (N1, N2, m) => {
    
    if (m < 0) return 0;
    else if (m > (N1+N2)) return 0;
    else {
        let sum = 0;

        for (let k = 0; k <= N2; ++k) {
            
            let fac1 = k;
            let fac2 = m - k;
            let fac3 = N2 - k;
            let fac4 = N1 - m + k;


            if (fac1 >= 0 && fac2 >= 0 && fac3 >= 0 && fac4 >= 0) {
                sum += Math.pow(-1, k) / (mathjs.factorial(fac1) * mathjs.factorial(fac2) * mathjs.factorial(fac3) * mathjs.factorial(fac4));   
            }
        }
        
        return mathjs.factorial(N1) * mathjs.factorial(N2) * sum;
    }
};

const isOdd = (x) => { return x & 1; };


const g0 = (alpha, beta, l1, l2, lamda) => {

    if (isOdd(l1 + alpha) || isOdd(l1 - alpha) || isOdd(l2 + beta) || isOdd(l2 - beta)) {
        return 0;
    } else {
        let sum = 0;

        for (let i = 0; i <= 2 * lamda; i = i+2) {
            if ((alpha + lamda -i) >= 0 && ((l1 - alpha) / 2 - lamda + (i / 2)) >= 0 && ((l1 + alpha) / 2 + lamda - (i / 2)) >= 0) {
                let counter = Math.pow(-1, i / 2) * mathjs.factorial(l1 + alpha + (2 * lamda) - i) * F(lamda, lamda, i);
                let divider = mathjs.factorial((l1 - alpha) / 2 - lamda + (i / 2)) * mathjs.factorial((l1 + alpha) / 2 + lamda - (i / 2)) * mathjs.factorial(alpha + lamda -i);

                sum += (counter / divider);
            }
        }

        let a = 1 / Math.pow(2, l1 + l2 + 1);
        let b = Math.sqrt(((2 * l1 + 1) * (mathjs.factorial(l1 - lamda)) * (2 * l2 + 1) * (mathjs.factorial(l2 - lamda))) / (mathjs.factorial(l1 + lamda) * mathjs.factorial(l2 + lamda)));
        let c = (mathjs.factorial(l1 + beta) * Math.pow(-1, (l1 - alpha) / 2 + (l2 - beta) / 2 - lamda)) / (mathjs.factorial((l2 - beta) / 2) * mathjs.factorial((l2 + beta) / 2) * mathjs.factorial(beta - lamda));

        return a * b * c * sum;
    }
} 

const pDot = (l1, l2, lamda, q, alpha, beta) => {
    return g0(alpha, beta, l1, l2, lamda) * F(alpha + lamda, beta - lamda, q);
};

const A = (ro, k = 0) => {
    
    let res;

    if (k == 0) {
        res = Math.exp(-1 * ro) / ro;
    } else if (k > 0) {
        res = (1 / ro) * ((k * A(ro, k - 1)) + Math.exp(-1 * ro));
    }

    return res;
};

const calcKTop = (kMax, ro, t, d = 16) => {

    let kTop;
    
    if (kMax == (ro * t)) {
        kTop =  Math.pow(kMax, 2);
    } else {
        let l = Math.log10(Math.abs(kMax / (ro * t)));
        kTop = Math.floor((d / Math.abs(l)) + kMax);
    }

    if (isOdd(kTop)) {
        kTop++;
    }

    return kTop;
};

const B = (ro, k, q) => {
    
    let res;
    const kMax = n1+n2+l1+l2;


    if (ro*t == 0) {
        res = (1 + Math.pow(-1, k)) / (k + 1);
    } else if (Math.abs(kMax / (ro * t)) < 1) {
        if (k == 0) {
            res = (1 / (ro *t)) * (Math.exp(ro * t) - Math.exp(-1 * ro * t));
        }
        else {
            res = (1 / (ro *t)) * (k * B(ro, k - 1, q) + Math.pow(-1, k) * Math.exp(ro * t) - Math.exp(-1 * ro * t));
        }
    } else {
        const kTop = calcKTop(kMax, ro, t);

        if (k == kTop) {
            res = (2 / (ro * t)) * Math.sinh(ro * t); 
        } else if (k < kTop) {
            
            const recB = B(ro, k + 1, q);



            res = (1 / (k + 1)) * ((ro * t * recB + Math.pow(-1, k) * Math.exp(ro * t) + Math.exp(-1 * ro * t)));
        } else {
            console.log('k > kTop!!!!!');
        }
    }

    return res;
};

const FAB = (N1, N2, q, ro, t) => {
    let sum = 0;

    for (let m = 0; m <= N1 + N2; ++m) {
        const currF = F(N1, N2, m);
        const currA = A(ro, N1 + N2 + q - m);
        const currB = B(ro, q + m, q);

        sum += currA * currB * currF;
    }

    return sum;
};

const A3 = (Xa) => {

    let loopsSum = 0;

    for (let alpha = -1*lamda; alpha <= l1; ++alpha) {
        for (let beta = lamda; beta <= l2; ++beta) {
            for (let q = 0; q <= alpha + beta; ++q) {
                const pd = pDot(l1, l2, lamda, q, alpha, beta);
                const currFAB = FAB(n1 - alpha, n2 - beta, q, ro(Xa), t); 
                loopsSum +=  pd * currFAB;
                
            }
        }
    }
    const rn = Math.pow(R(Xa), N);

    console.log(`loopsSum = ${loopsSum}, Htag = ${Htag()}, R^N = ${rn}`);

    let A3 = Htag() * rn * loopsSum;

    return A3;
};

// ====================================================================

const manualDerivativeOneParam = (f, a, method, h = 0.001) => {

    if (method == 'central') {
        let first = f(a+h);
        console.log(`first = ${first}`);
        let second = f(a-h);
        console.log(`second = ${second}`);
        return (first - second) / (2*h);
    } else if (method == 'forward') {
        return (f(a + h) - f(a)) / h;
    } else if (method == 'backward') {
        return (f(a) - f(a - h)) / h;
    } else throw 'Method must be "central", "forward" or "backward"';
};

const manualDerivativeTwoParams = (f, a, b, method, h = 0.00000000001) => {

    if (method == 'central') {
        return (f(a + h, b) - f(a - h, b)) / (2*h);
    } else if (method == 'forward') {
        return (f(a + h, b) - f(a, b)) / h;
    } else if (method == 'backward') {
        return (f(a, b) - f(a - h, b)) / h;
    } else throw 'Method must be "central", "forward" or "backward"';
};

const manualDerivativeFourParams = (f, a, b, c, d, method, h = 0.00000000001) => {

    if (method == 'central') {
        return (f(a + h, b, c, d) - f(a - h, b, c, d)) / (2*h);
    } else if (method == 'forward') {
        return (f(a + h, b, c, d) - f(a, b, c, d)) / h;
    } else if (method == 'backward') {
        return (f(a, b, c, d) - f(a - h, b, c, d)) / h;
    } else throw 'Method must be "central", "forward" or "backward"';
};

const Htag = () => {
    return H() * Math.pow(h1/2, N);
};

const resultToCompare = (Xa) => {
    
    let loopsSum1 = 0;

    for (let alpha = -1*lamda; alpha <= l1; ++alpha) {
        for (let beta = lamda; beta <= l2; ++beta) {
            for (let q = 0; q <= alpha + beta; ++q) {
                
                const currpDot = pDot(l1, l2, lamda, q, alpha, beta);
                const currFAB = FAB(n1 - alpha, n2 - beta, q, ro(Xa), t);

                //console.log(`pDot = ${currpDot}, FAB = ${currFAB}`);
                loopsSum1 += currpDot * currFAB;
            }
        }
    }

    let loopsSum2 = 0;

    for (let alpha = -1*lamda; alpha <= l1; ++alpha) {
        for (let beta = lamda; beta <= l2; ++beta) {
            for (let q = 0; q <= alpha + beta; ++q) {

                let Q2 = 0;
                let N1 = n1 - alpha;
                let N2 = n2 - beta;

                for (let m = 0; m <= N1 + N2; ++m) {
                    Q2 += F(N1, N2, m) * ((((Xa-Xb)/R(Xa))* -1*h1/2 *
                        A(ro(Xa), N1 + N2 + q - m + 1) * B(ro(Xa), q + m, q)) + 
                        ((((Xa-Xb)/R(Xa))* -1*h1/2) * t * A(ro(Xa), N1 + N2 + q - m) * B(ro(Xa), q + m + 1, q)));
                }


                loopsSum2 += pDot(l1, l2, lamda, q, alpha, beta) * Q2;
            }
        }
    }

    return (Htag() * N * (Xa-Xb) * Math.pow(R(Xa),N-2) * loopsSum1) + 
        (Htag() * Math.pow(R(Xa), N) * loopsSum2);
};

// ==================================

const main = (Xa, h = 0.01) => {
    console.log(`A3 = ${A3(Xa)}`);
    const A3derivative = manualDerivativeOneParam(A3, Xa, 'central', h);
    const compare = resultToCompare(Xa);

    console.log(`A3 derivative = ${A3derivative}`);
    console.log(`result to compare = ${compare}`);
};


module.exports = main;