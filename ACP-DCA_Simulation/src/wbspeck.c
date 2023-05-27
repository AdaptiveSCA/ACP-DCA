#include "wbspeck.h"

void Key_expand32(u16 key[4], u16 roundkey[22])
{
    u16 i, b = key[3];
    u16 a[SPECK_KEY_LEN - 1];

    for(i = 0; i < (SPECK_KEY_LEN - 1); i++)
    {
        a[i] = key[i + 1];
    }
    roundkey[0] = b;  
    for(i = 0; i < SPECK_ROUNDS - 1; i++)
    {
        R32(a[i % (SPECK_KEY_LEN - 1)], b, i);
        roundkey[i + 1] = b;
    }
}
void Key_expand64(u32 key[4], u32 roundkey[22])
{
    u32 i, b = key[3];
    u32 a[SPECK_KEY_LEN - 1];

    for(i = 0; i < (SPECK_KEY_LEN - 1); i++)
    {
        a[i] = key[2 - i];
    }
    roundkey[0] = b;  
    for(i = 0; i < SPECK_ROUNDS - 1; i++)
    {
        R64(a[i % (SPECK_KEY_LEN - 1)], b, i);
        roundkey[i + 1] = b;
    }
}
void Key_expand128(u64 key[4], u64 roundkey[22])
{
    u64 i, b = key[3];
    u64 a[SPECK_KEY_LEN - 1];

    for(i = 0; i < (SPECK_KEY_LEN - 1); i++)
    {
        a[i] = key[2 - i];
    }
    roundkey[0] = b;  
    for(i = 0; i < SPECK_ROUNDS - 1; i++)
    {
        R128(a[i % (SPECK_KEY_LEN - 1)], b, i);
        roundkey[i + 1] = b;
    }
}

u16 ROR32(u16 x, int r)
{
    return (x >> r) | (x << ((sizeof(u16) * 8) - r));
}
u16 ROL32(u16 x, int r)
{
    return (x << r) | (x >> ((sizeof(u16) * 8) - r));
}
u32 ROR64(u32 x, int r)
{
    return (x >> r) | (x << ((sizeof(u32) * 8) - r));
}
u32 ROL64(u32 x, int r)
{
    return (x << r) | (x >> ((sizeof(u32) * 8) - r));
}
u64 ROR128(u64 x, int r)
{
    return (x >> r) | (x << ((sizeof(u64) * 8) - r));
}
u64 ROL128(u64 x, int r)
{
    return (x << r) | (x >> ((sizeof(u64) * 8) - r));
}

void ACP_DCA32()
{
    u16 x0, x1, y0, y1, z0, z1;
    int k, x, j, s, ts, l, b;
    u32 state;
    u16 key[4] = {0xe318, 0x2610, 0x7a08, 0x1918};
    u16 roundkey[SPECK_ROUNDS];
    Key_expand32(key, roundkey);
    ////
    int f = 0, g = 0;
    static const int ktb16 = 0x0000;
    static const int  kte16 = 0xffff;
    int key_count = 0;

    Aff32 Aff, Aff_inv;
    genaffinepairM32(&Aff, &Aff_inv);
    u32 map[pt]; 
    
    u16 key_guess[10000] = {0};
    int k_count = 0;
    double k_max = 0.0;
    int knum = 0;
    for(k = ktb16; k < kte16; k++) // key
    {
        for(x = 0; x < pt; x++) // adaptive inputs
        {
            x0 = x;
            x1 = 0x00ff;
            
            y1 = ROR32(ROL32(x0 - x1, 7) ^ x1, 2);
            y0 = ROL32((ROL32(x0 - x1, 7) ^ k) - y1, 7);
            
            z1 = ROL32(y1, 2) ^ (ROR32(y0, 7) + y1) ^ roundkey[0];
            z0 = ROR32((ROR32(y0, 7) + y1) ^ roundkey[0], 7) + z1;

            state = (z0 << 16) | z1;
            map[x] = state;
            map[x] = affineU32(Aff, state); // affine encoding
        }
        double k_score = 0.0;
        double L_score[pt] = {0.0};
        int l_count = 0;
        double l_max = 0.0;
        for(l = 1; l < pt; l++)
        {  
            double score = 0.0;
            double j_max = 0.0;
            for(j = 0; j < tb32; j++) // samples of traces
            {
                int Nf0 = 0, Nf1 = 0, Ng0 = 0, Ng1 = 0, N00 = 0, N01 = 0, N10 = 0, N11 = 0;
                for(x = 0; x < pt; x++)
                { 
                    if(map[x] & idM32[j]) ts = 1; // traces
                    else ts = 0;
                    ////// correlation
                    f = ts;
                    if(f) Nf1++;
                    else Nf0++;
                    
                    g = xor[l & x];
                    if(g) Ng1++;
                    else Ng0++;
                    
                    if(f == 1 && g == 1) N11++;
                    else if(f == 0 & g == 0) N00++;
                    else if(f == 1 & g == 0) N10++;
                    else N01++;
                }
                if(Nf1 && Nf0 && Ng1 && Ng0) score = abs((N11 * N00 - N10 * N01)) * 1.0 / (sqrt(Nf1) * sqrt(Nf0) * sqrt(Ng1) * sqrt(Ng0));
                else score = 0.0;
                if(score > j_max) j_max = score;
            }
            L_score[l] = j_max;
            if(L_score[l] > l_max) l_max = L_score[l];
        }
        k_score = l_max;
        if(fabs(k_score - k_max) < EPS)
        {
            for(l = 1; l < pt; l++)
            {
                if(fabs(L_score[l] - l_max) < EPS) // each encoding: key guess
                {
                    l_count++;
                }
            }
            if(l_count > k_count)
            {
                k_count = l_count;
                knum = 0;
                key_guess[knum] = k;
                knum++;
                printf("key guess: %.2x, encoding count: %d, correlation: %f\n", k, l_count, k_score);
            }
            else if(l_count == k_count)
            {
                key_guess[knum] = k;
                knum++;
                printf("key guess: %.2x, encoding count: %d, correlation: %f\n", k, l_count, k_score);
            }
        }
        else if(k_score > k_max) 
        {
            k_max = k_score; 
            knum = 0;
            for(l = 1; l < pt; l++)
            {
                if(fabs(L_score[l] - l_max) < EPS) // each encoding: key guess
                {
                    l_count++;
                }
            }
            k_count = l_count;
            key_guess[knum] = k;
            knum++;
            printf("key guess: %.2x, encoding count: %d, correlation: %f\n", k, l_count, k_score);
        }
    }
    printf("------\n");  
    for(k = 0; k < knum; k++)
    {
        printf("simulation SPECK32 recoverd key: %.2x, correlation: %f, encoding count: %d\n", key_guess[k], k_max, k_count);
        key_count++;
    }
    printf("simulation SPECK32 recoverd key count: %d\n", key_count);
}

void ACP_DCA64()
{
    u32 x0, x1, y0, y1, z0, z1;
    int k, x, j, s, ts, l, b;
    u64 state;
    u32 key[4] = {0xe318e318, 0x26102610, 0x7a087a08, 0x13121110};
    u32 roundkey[SPECK_ROUNDS];
    Key_expand64(key, roundkey);
    ////
    int f = 0, g = 0;
    static const int ktb32 = 0x13100000;
    static const int  kte32 = 0x131fffff;
    int key_count = 0;

    Aff64 Aff, Aff_inv;
    genaffinepairM64(&Aff, &Aff_inv);
    u64 map[pt]; 
    
    u32 key_guess[10000] = {0};
    int k_count = 0;
    double k_max = 0.0;
    int knum = 0;
    for(k = ktb32; k < kte32; k++) // key
    {
        for(x = 0; x < pt; x++) // adaptive inputs
        {
            x0 = x;
            x1 = 0x00ff00ff;
            
            y1 = ROR64(ROL64(x0 - x1, 8) ^ x1, 3);
            y0 = ROL64((ROL64(x0 - x1, 8) ^ k) - y1, 8);
            
            z1 = ROL64(y1, 3) ^ (ROR64(y0, 8) + y1) ^ roundkey[0];
            z0 = ROR64((ROR64(y0, 8) + y1) ^ roundkey[0], 8) + z1;

            state = ((u64)(z0) << 32) | z1;
            map[x] = state;
            map[x] = affineU64(Aff, state); // affine encoding
        }
        double k_score = 0.0;
        double L_score[pt] = {0.0};
        int l_count = 0;
        double l_max = 0.0;
        for(l = 1; l < pt; l++)
        {  
            double score = 0.0;
            double j_max = 0.0;
            for(j = 0; j < tb64; j++) // samples of traces
            {
                int Nf0 = 0, Nf1 = 0, Ng0 = 0, Ng1 = 0, N00 = 0, N01 = 0, N10 = 0, N11 = 0;
                for(x = 0; x < pt; x++)
                { 
                    if(map[x] & idM64[j]) ts = 1; // traces
                    else ts = 0;
                    ////// correlation
                    f = ts;
                    if(f) Nf1++;
                    else Nf0++;
                    
                    g = xor[l & x];
                    if(g) Ng1++;
                    else Ng0++;
                    
                    if(f == 1 && g == 1) N11++;
                    else if(f == 0 & g == 0) N00++;
                    else if(f == 1 & g == 0) N10++;
                    else N01++;
                }
                if(Nf1 && Nf0 && Ng1 && Ng0) score = abs((N11 * N00 - N10 * N01)) * 1.0 / (sqrt(Nf1) * sqrt(Nf0) * sqrt(Ng1) * sqrt(Ng0));
                else score = 0.0;
                if(score > j_max) j_max = score;
            }
            L_score[l] = j_max;
            if(L_score[l] > l_max) l_max = L_score[l];
        }
        k_score = l_max;
        if(fabs(k_score - k_max) < EPS)
        {
            for(l = 1; l < pt; l++)
            {
                if(fabs(L_score[l] - l_max) < EPS) // each encoding: key guess
                {
                    l_count++;
                }
            }
            if(l_count > k_count)
            {
                k_count = l_count;
                knum = 0;
                key_guess[knum] = k;
                knum++;
                printf("key guess: %.2x, encoding count: %d, correlation: %f\n", k, l_count, k_score);
            }
            else if(l_count == k_count)
            {
                key_guess[knum] = k;
                knum++;
                printf("key guess: %.2x, encoding count: %d, correlation: %f\n", k, l_count, k_score);
            }
        }
        else if(k_score > k_max) 
        {
            k_max = k_score; 
            knum = 0;
            for(l = 1; l < pt; l++)
            {
                if(fabs(L_score[l] - l_max) < EPS) // each encoding: key guess
                {
                    l_count++;
                }
            }
            k_count = l_count;
            key_guess[knum] = k;
            knum++;
            printf("key guess: %.2x, encoding count: %d, correlation: %f\n", k, l_count, k_score);
        }
    }
    printf("------\n");  
    for(k = 0; k < knum; k++)
    {
        printf("simulation SPECK64 recoverd key: %.2x, correlation: %f, encoding count: %d\n", key_guess[k], k_max, k_count);
        key_count++;
    }
    printf("simulation SPECK64 recoverd key count: %d\n", key_count);
}

void ACP_DCA128()
{
    u64 x0, x1, y0, y1, z0, z1;
    int x, j, s, ts, l, b;
    u64 k;
    V128 state1;
    V128 state2;
    u64 key[4] = {0xe318e318, 0x26102610, 0x7a087a08, 0x0f0e0d0c0b0a0908};
    u64 roundkey[SPECK_ROUNDS];
    Key_expand128(key, roundkey);
    ////
    int f = 0, g = 0;
    static const u64 ktb64 = 0x0f0e0d0c0b0a0000;
    static const u64 kte64 = 0x0f0e0d0c0b0fffff;
    int key_count = 0;

    Aff128 Aff, Aff_inv;
    genaffinepairM128(&Aff, &Aff_inv);
    u64 map[pt][2]; 
    
    u64 key_guess[10000] = {0};
    int k_count = 0;
    double k_max = 0.0;
    int knum = 0;
    for(k = ktb64; k < kte64; k++) // key
    {
        for(x = 0; x < pt; x++) // adaptive inputs
        {
            x0 = x;
            x1 = 0x00ff00ff00ff00ff;
            
            y1 = ROR128(ROL128(x0 - x1, 8) ^ x1, 3);
            y0 = ROL128((ROL128(x0 - x1, 8) ^ k) - y1, 8);
            
            z1 = ROL128(y1, 3) ^ (ROR128(y0, 8) + y1) ^ roundkey[0];
            z0 = ROR128((ROR128(y0, 8) + y1) ^ roundkey[0], 8) + z1;

            state1.V[0] = z0; 
            state1.V[1] = z1;
            MatMulVecM128(Aff.Mat, state1, &state2);
            map[x][0] = state2.V[0] ^ Aff.Vec.V[0];//add
            map[x][1] = state2.V[1] ^ Aff.Vec.V[1];
        }
        double k_score = 0.0;
        double L_score[pt] = {0.0};
        int l_count = 0;
        double l_max = 0.0;
        for(l = 1; l < pt; l++)
        {  
            double score = 0.0;
            double j_max = 0.0;
            for(j = 0; j < tb128; j++) // samples of traces
            {
                int Nf0 = 0, Nf1 = 0, Ng0 = 0, Ng1 = 0, N00 = 0, N01 = 0, N10 = 0, N11 = 0;
                for(x = 0; x < pt; x++)
                { 
                    if(j < 64)
                    {    
                        if(map[x][0] & idM64[j]) ts = 1; // traces
                        else ts = 0;
                    }
                    else
                    {    
                        if(map[x][1] & idM64[j - 64]) ts = 1; // traces
                        else ts = 0;
                    }
                    ////// correlation
                    f = ts;
                    if(f) Nf1++;
                    else Nf0++;
                    
                    g = xor[l & x];
                    if(g) Ng1++;
                    else Ng0++;
                    
                    if(f == 1 && g == 1) N11++;
                    else if(f == 0 & g == 0) N00++;
                    else if(f == 1 & g == 0) N10++;
                    else N01++;
                }
                if(Nf1 && Nf0 && Ng1 && Ng0) score = abs((N11 * N00 - N10 * N01)) * 1.0 / (sqrt(Nf1) * sqrt(Nf0) * sqrt(Ng1) * sqrt(Ng0));
                else score = 0.0;
                if(score > j_max) j_max = score;
            }
            L_score[l] = j_max;
            if(L_score[l] > l_max) l_max = L_score[l];
        }
        k_score = l_max;
        if(fabs(k_score - k_max) < EPS)
        {
            for(l = 1; l < pt; l++)
            {
                if(fabs(L_score[l] - l_max) < EPS) // each encoding: key guess
                {
                    l_count++;
                }
            }
            if(l_count > k_count)
            {
                k_count = l_count;
                knum = 0;
                key_guess[knum] = k;
                knum++;
                printf("key guess: %llx, encoding count: %d, correlation: %f\n", k, l_count, k_score);
            }
            else if(l_count == k_count)
            {
                key_guess[knum] = k;
                knum++;
                printf("key guess: %llx, encoding count: %d, correlation: %f\n", k, l_count, k_score);
            }
        }
        else if(k_score > k_max) 
        {
            k_max = k_score; 
            knum = 0;
            for(l = 1; l < pt; l++)
            {
                if(fabs(L_score[l] - l_max) < EPS) // each encoding: key guess
                {
                    l_count++;
                }
            }
            k_count = l_count;
            key_guess[knum] = k;
            knum++;
            printf("key guess: %llx, encoding count: %d, correlation: %f\n", k, l_count, k_score);
        }
    } 
    printf("------\n");  
    for(k = 0; k < knum; k++)
    {
        printf("simulation SPECK128 recoverd key: %llx, correlation: %f, encoding count: %d\n", key_guess[k], k_max, k_count);
        key_count++;
    }
    printf("simulation SPECK128 recoverd key count: %d\n", key_count);
}