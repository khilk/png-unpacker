#include <iostream>
#include <array>
#include <vector>
#include <cstdio>
#include <future>
#include <chrono>

using namespace std;

void GetValue(vector<unsigned>& bit_data, int& pointer, int& L, int& B, unsigned& value){
    if (((value << 1u) + 1) <= 191 && value << 1u >= 48) { // первый случай, длина 8 бит
        value = (value << 1u) + bit_data[pointer++];
        L = 0;
        B = 48;
    } else if ((value << 2u) + 3 <= 511 && value << 2u >= 400) { // второй случай, длина 9 бит
        value = (value << 2u) + (bit_data[pointer] << 1u) + bit_data[pointer + 1];
        pointer += 2;
        L = 144;
        B = 400;
    } else if (value >= 0 && value <= 23) { // третий случай, длина 7 бит, смотрим таблицы длин и смещений
        L = 256;
        B = 0;
    } else if (((value << 1u) + 1) <= 199 &&
               value << 1u >= 192) { // четвёртый случай, длина 8 бит, смотрим таблицу
        value = (value << 1u) + bit_data[pointer++];
        L = 192;
        B = 280;
    } else {
        throw runtime_error( "Something wrong. Terminate");
    }
}

void GetLength(vector<unsigned>& bit_data, int& pointer, unsigned& value_length, unsigned V){
    if(V <= 264){
        value_length = V - 254;
    }
    else{
        unsigned add_bits = (V - 261) / 4;
        unsigned constanta = 10;
        for(unsigned i = 1; i < add_bits; ++i){
            constanta += 4u * (1u << i);
        }
        constanta += ((V - 265) % 4) * (2u << (add_bits - 1)) + 1;
        unsigned add = 0;
        for(int i = 0; i < add_bits; ++i){
            add += bit_data[pointer + i] << unsigned (i);
        }
        pointer += add_bits;
        value_length = constanta + add;
    }
}

void GetShift(vector<unsigned>& bit_data, int& pointer, unsigned& shift_value, unsigned shift){
    if(shift <= 3){
        shift_value = shift + 1;
    }
    else{
        unsigned add_bits = (shift - 2) / 2;
        unsigned constanta = 4;
        for(unsigned i = 1; i < add_bits; ++i){
            constanta += 2u * (1u << i);
        }
        if(shift % 2 == 1){
            constanta += 2u << (add_bits - 1);
        }
        constanta += 1;
        unsigned add = 0;
        for(int i = 0; i < add_bits; ++i){
            add += bit_data[pointer + i] << unsigned (i);
        }
        pointer += add_bits;
        shift_value = constanta + add;
    }
}

unsigned Adler32(vector<unsigned>& decompressed_data){
    unsigned A = 1, B = 0;
    int n = decompressed_data.size();
    for(int i = 0; i < n; ++i){
        A += decompressed_data[i];
        B += (n - i) * decompressed_data[i] + 1;
    }
    A %= 65521;
    B %= 65521;
    return B * 65536 + A;
}

void Decode(vector<unsigned>& bit_data, vector<unsigned>& decompressed_data, int w, int h){
    int pointer = 0 , bfinal = 0;
    while(bfinal != 1) {
        bfinal = bit_data[pointer];
        pointer++;
        unsigned type = (bit_data[pointer] << 1u) + bit_data[pointer + 1];
        if (type == 1) { // динамические коды Хаффмана
            throw runtime_error("This program deals only with fixed Huffman codes! Terminate");
        }
        pointer += 2;
        if(type == 2) { // фиксированные коды Хаффмана
            // смотрим таблицу фиксированных кодов Хаффмана
            int L = 0, B = 0;
            unsigned value = 0;
            unsigned j = 0;
            while (true) {
                value = 0;
                for (int i = 0; i < 7; ++i) {
                    value += bit_data[pointer + i] << unsigned(6 - i);
                }
                pointer += 7;
                GetValue(bit_data, pointer, L, B, value);
                unsigned V = value - B + L;
                if (V == 256) {
                    break;
                }
                if (0 <= V && V <= 255) {
                    decompressed_data.push_back(V);
                    j++;
                } else {
                    unsigned value_length = 0, shift = 0, shift_value = 0;
                    GetLength(bit_data, pointer, value_length, V);
                    for (int i = 0; i < 5; ++i) {
                        shift += bit_data[pointer + i] << unsigned(4 - i);
                    }
                    pointer += 5;
                    GetShift(bit_data, pointer, shift_value, shift);
                    int k = 0;
                    unsigned index = decompressed_data.size() - shift_value;
                    for (int i = 0; i < value_length; ++i) {
                        auto new_value = decompressed_data[index + k];
                        decompressed_data.push_back(new_value);
                        k++;
                        if (k == shift_value) {
                            k = 0;
                        }
                    }
                    j += value_length;
                }
            }
        }
        if(type == 0){ // без сжатия
            while(pointer % 8 != 0){ // выравнивание на границу байта
                pointer++;
            }
            int LEN = 0;
            for (int i = 0; i < 16; ++i) {
                LEN += bit_data[pointer + i] << unsigned(i);
            }
            pointer += 32;
            for(int i = 0; i < LEN; ++i){
                int value = 0;
                for(int j = 0; j < 8; ++j){
                    value += bit_data[pointer + j] << unsigned(j);
                }
                pointer += 8;
                decompressed_data.push_back(value);
            }
        }
    }
    unsigned control = 0, n = bit_data.size() - 32; // проверка контрольных сумм
    int k = 0;
    for(int j = 0; j < 4; ++j){
        for(int i = 0; i < 8; ++i){
        control += bit_data[n + 7 - i] << unsigned(31 - k);
        k++;
        }
        n += 8;
    }
    if(control != Adler32(decompressed_data)) {
        throw runtime_error("Control sum incorrect");
    }
}

void Inflate(FILE *source, std::vector<unsigned>& decompressed_data, int w, int h)
{
    fseek(source, 33, SEEK_SET); //skip file header and IHDR
    uint8_t end_id[4] = {0x49, 0x45, 0x4E, 0x44};
    uint8_t chunk_id[4] = {0x49, 0x44, 0x41, 0x54};
    bool is_end = false;
    while(!is_end) {
        int _length;
        while (true) {// search for IDAT chunk
            _length = 0;
            uint8_t length[4];
            int read = fread(length, sizeof(uint8_t), 4, source);
            if (read != 4) {
                throw runtime_error("Error reading file! Terminate");
            }
            for (int i = 0; i < 4; i++) {
                _length += length[i] << unsigned (8 * (3 - i));
            }
            uint8_t b;
            int incorrect = 0;
            is_end = true;
            for (int i = 0; i < 4; i++) {
                fread(&b, sizeof(uint8_t), 1, source);
                if (b != chunk_id[i]) {
                    incorrect++;
                }
                if (b != end_id[i]){
                    is_end = false;
                }
            }
            if (incorrect == 0 || is_end) {
                break;
            }
            fseek(source, _length + 4, SEEK_CUR);
        }
        if(!is_end && _length != 0) { // if not end -> decompress
            // filter - none
            uint8_t second;
            fseek(source, 1, SEEK_CUR);
            fread(&second, sizeof(uint8_t), 1, source);
            int k = 0;
            while (second != 0 && k != 6) {
                int bit_value = int(second) % 2;
                second = int(second) / 2;
                k++;
                if(k == 6 && bit_value == 1){
                    throw runtime_error("this program doesn't deal with dictionaries");
                }
            }
            vector<unsigned> bit_data;
            bit_data.reserve(_length * 8);
            for(int j = 0; j < _length - 2; ++j) {
                uint8_t b;
                fread(&b, sizeof(uint8_t), 1, source);
                int i = 0;
                while (b != 0) {
                    bit_data.push_back(int(b) % 2);
                    b = int(b) / 2;
                    i++;
                }
                while (i != 8) {
                    bit_data.push_back(0);
                    i++;
                }
            }
            Decode(bit_data, decompressed_data, w, h);
            fseek(source, 4, SEEK_CUR);
        }
    }
}

void CheckFormat(FILE *source){
    fseek(source, 1, SEEK_SET);
    vector<uint8_t> png_id = { 0x50, 0x4E, 0x47};
    for(auto& item : png_id) {
        uint8_t n;
        fread(&n, sizeof(uint8_t), 1, source);
        if (n != item) {
            throw runtime_error("Wrong format! Terminate");
        }
    }
    fseek(source, 24, SEEK_SET);
    uint8_t bit, color, compression, filter, interlace;
    fread(&bit, sizeof(uint8_t), 1, source);
    fread(&color, sizeof(uint8_t), 1, source);
    fread(&compression, sizeof(uint8_t), 1, source);
    fread(&filter, sizeof(uint8_t), 1, source);
    fread(&interlace, sizeof(uint8_t), 1, source);
    if(filter != 0){ // filter is always equal to 0 in png
        throw runtime_error("filter incorrect");
    }
    if(interlace != 0){ // this program deals only with images with no interlace
        throw runtime_error("interlace method  incorrect");
    }
    if(compression != 0){ // always 0, Deflate
        throw runtime_error("compression method  incorrect");
    }
}

void Gistogramma(vector<unsigned>& decompressed_data, array<int, 256>& gistogramma, uint32_t w, uint32_t h, int k){ // for 8 bit rgb
    vector<future<array<int, 256>>> futures;
    //int k = 1; // number of parts
    int last = h;
    int begin = 0, end = 0;
    auto start = chrono::high_resolution_clock::now();
    for(int i = 0; i < k; ++i){
        int part_length = h % k == 0 ? h / k : (h / k + 1);
        int mult = part_length < last ? part_length : last;
        end = begin + (3 * int(w) + 1) * mult;
        futures.push_back(
                async([&decompressed_data, begin, end, w]{
                    array<int, 256> g;
                    g.fill(0);
                    int i = begin + 1, l = 0 , j = 1;
                    while(i < end){
                        j++;
                        int value = int(0.3 * decompressed_data[i]
                                + 0.59 * decompressed_data[i + 1] + 0.11 * decompressed_data[i + 2]);
                        i += 3;
                        l++;
                        if(l == w){
                            i++;
                            l = 0;
                        }
                        g[value]++;
                    }
                return g;})
                );
        begin = end;
        last -= part_length;
    }
    for(auto& item : futures){
        int i = 0;
        for(auto& v : item.get()){
            gistogramma[i] += v;
            i++;
        }
    }
    auto finish = chrono::high_resolution_clock::now();
    chrono::duration<double> elapsed = finish - start;
    cout << "Elapsed time: " << elapsed.count() << endl;
}

int main() {
    FILE* source = fopen("/home/daria/Downloads/example3.png", "rb");
    if(source == nullptr) {
        cout << "Error opening file! Terminate" << endl;
        return 1;
    }
    FILE* to_write = fopen("/home/daria/Downloads/text.txt", "wr");
    if(to_write == nullptr) {
        cout << "Error opening file! Terminate" << endl;
        return 1;
    }
    try {
        CheckFormat(source);
        fseek(source, 16, SEEK_SET);
        uint32_t h = 0, w = 0;
        array<uint8_t, 4> width{0};
        array<uint8_t, 4> height{0};
        for (auto &item : width) {
            fread(&item, sizeof(uint8_t), 1, source);
        }
        for (auto &item : height) {
            fread(&item, sizeof(uint8_t), 1, source);
        }
        for (int i = 0; i < 4; i++) {
            w += width[i] << unsigned(8 * (3 - i));
        }
        for (int i = 0; i < 4; i++) {
            h += height[i] << unsigned(8 * (3 - i));
        }
        vector<unsigned> decompressed_data;
        decompressed_data.reserve(w * h * 6 + h);
        Inflate(source, decompressed_data, w, h);
        cout << "decompressed data:" << endl;
        for(auto& item : decompressed_data){
            cout << item << " ";
        }
        cout << endl << "size: " << decompressed_data.size() << endl;
        array<int, 256> gistogramma;
        gistogramma.fill(0);
        /*vector<unsigned > test(69000000);
        for(int i = 0; i < test.size(); ++i){
            test[i] = i % 256;
        }*/
        Gistogramma(decompressed_data, gistogramma, w, h, 1);
        int sum = 0;
        for(auto& item : gistogramma){
            cout << item << " ";
            //fprintf(to_write, "%d",  item);
            //fprintf(to_write, "%s",  " ");
            //sum += item;
        }
        gistogramma.fill(0);
        cout << endl;
        Gistogramma(decompressed_data, gistogramma, w, h, 4);
        for(auto& item : gistogramma){
            cout << item << " ";
        }
        //cout << "sum: " << sum << endl;
    }
    catch (exception& ex) {
        cout << ex.what() << endl;
    }
    return 0;
}
