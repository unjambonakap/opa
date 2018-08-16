#include <opa_common.h>
#include <opa/crypto/misc.h>

using namespace opa::crypto;
using namespace std;

int main() {
    int n = 3;
    u32 tb[] = { 6, 24, 1, 13, 16, 10, 20, 17, 15 };

    vector<u32> m;
    m.assign(tb, tb + n * n);
    HillCipher *x = HillCipher::get(26, n);
    x->setKey(m);

    u32 input_tb[] = { 0, 2, 19 };
    vector<u32> input;
    input.assign(input_tb, input_tb + n);

    vector<u32> cipher = x->encrypt(input);
    out(cipher);

    vector<u32> res = x->decrypt(cipher);
    out(res);
    out(input);

    return 0;
}
