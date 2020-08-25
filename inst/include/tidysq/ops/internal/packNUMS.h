#ifndef TIDYSQ_PACKNUMS_H
#define TIDYSQ_PACKNUMS_H

#include <stdexcept>

#include "../../types/general.h"
#include "util.h"

namespace tidysq::internal {
    template<InternalType INTERNAL_IN,
            ProtoType PROTO_IN,
            InternalType INTERNAL_OUT>
    Sequence<INTERNAL_OUT> packNUMS2(const SequenceProto<INTERNAL_IN, PROTO_IN> &unpacked,
                                     const Alphabet &alphabet) {
        lensq outByte = 0;
        lensq i = 0;
        Sequence<INTERNAL_OUT> packed = reserveSpaceForPacked<INTERNAL_IN, PROTO_IN, INTERNAL_OUT>(unpacked, alphabet);
        for (; i + 8 <= unpacked.size(); i += 8) {
            packed[outByte] = (unpacked[i]) |
                              (unpacked[i + 1] << 2) |
                              (unpacked[i + 2] << 4) |
                              (unpacked[i + 3] << 6);
            packed[outByte + 1] = (unpacked[i + 4]) |
                                  (unpacked[i + 5] << 2) |
                                  (unpacked[i + 6] << 4) |
                                  (unpacked[i + 7] << 6);
            outByte += 2;
        }
        switch (unpacked.size() - i) {
            case 7:
                packed[outByte] = (unpacked[i]) |
                                  (unpacked[i + 1] << 2) |
                                  (unpacked[i + 2] << 4) |
                                  (unpacked[i + 3] << 6);
                packed[outByte + 1] = (unpacked[i + 4]) |
                                      (unpacked[i + 5] << 2) |
                                      (unpacked[i + 6] << 4);
                break;
            case 6:
                packed[outByte] = (unpacked[i]) |
                                  (unpacked[i + 1] << 2) |
                                  (unpacked[i + 2] << 4) |
                                  (unpacked[i + 3] << 6);
                packed[outByte + 1] = (unpacked[i + 4]) |
                                      (unpacked[i + 5] << 2);
                break;
            case 5:
                packed[outByte] = (unpacked[i]) |
                                  (unpacked[i + 1] << 2) |
                                  (unpacked[i + 2] << 4) |
                                  (unpacked[i + 3] << 6);
                packed[outByte + 1] = (unpacked[i + 4]);
                break;
            case 4:
                packed[outByte] = (unpacked[i]) |
                                  (unpacked[i + 1] << 2) |
                                  (unpacked[i + 2] << 4) |
                                  (unpacked[i + 3] << 6);
                break;
            case 3:
                packed[outByte] = (unpacked[i]) |
                                  (unpacked[i + 1] << 2) |
                                  (unpacked[i + 2] << 4);
                break;
            case 2:
                packed[outByte] = (unpacked[i]) |
                                  (unpacked[i + 1] << 2);
                break;
            case 1:
                packed[outByte] = (unpacked[i]);
                break;
        }
        return packed;
    }

    template<InternalType INTERNAL_IN,
            ProtoType PROTO_IN,
            InternalType INTERNAL_OUT>
    Sequence<INTERNAL_OUT> packNUMS3(const SequenceProto<INTERNAL_IN, PROTO_IN> &unpacked,
                                     const Alphabet &alphabet) {
        lensq outByte = 0;
        lensq i = 0;
        Sequence<INTERNAL_OUT> packed = reserveSpaceForPacked<INTERNAL_IN, PROTO_IN, INTERNAL_OUT>(unpacked, alphabet);
        for (; i + 8 <= unpacked.size(); i += 8) {
            packed[outByte] = (unpacked[i]) |
                              (unpacked[i + 1] << 3) |
                              (unpacked[i + 2] << 6);
            packed[outByte + 1] = (unpacked[i + 2] >> 2) |
                                  (unpacked[i + 3] << 1) |
                                  (unpacked[i + 4] << 4) |
                                  (unpacked[i + 5] << 7);
            packed[outByte + 2] = (unpacked[i + 5] >> 1) |
                                  (unpacked[i + 6] << 2) |
                                  (unpacked[i + 7] << 5);
            outByte += 3;
        }
        switch (unpacked.size() - i) {
            case 7:
                packed[outByte] = (unpacked[i]) |
                                  (unpacked[i + 1] << 3) |
                                  (unpacked[i + 2] << 6);
                packed[outByte + 1] = (unpacked[i + 2] >> 2) |
                                      (unpacked[i + 3] << 1) |
                                      (unpacked[i + 4] << 4) |
                                      (unpacked[i + 5] << 7);
                packed[outByte + 2] = (unpacked[i + 5] >> 1) |
                                      (unpacked[i + 6] << 2);
                break;
            case 6:
                packed[outByte] = (unpacked[i]) |
                                  (unpacked[i + 1] << 3) |
                                  (unpacked[i + 2] << 6);
                packed[outByte + 1] = (unpacked[i + 2] >> 2) |
                                      (unpacked[i + 3] << 1) |
                                      (unpacked[i + 4] << 4) |
                                      (unpacked[i + 5] << 7);
                packed[outByte + 2] = (unpacked[i + 5] >> 1);
                break;
            case 5:
                packed[outByte] = (unpacked[i]) |
                                  (unpacked[i + 1] << 3) |
                                  (unpacked[i + 2] << 6);
                packed[outByte + 1] = (unpacked[i + 2] >> 2) |
                                      (unpacked[i + 3] << 1) |
                                      (unpacked[i + 4] << 4);
                break;
            case 4:
                packed[outByte] = (unpacked[i]) |
                                  (unpacked[i + 1] << 3) |
                                  (unpacked[i + 2] << 6);
                packed[outByte + 1] = (unpacked[i + 2] >> 2) |
                                      (unpacked[i + 3] << 1);
                break;
            case 3:
                packed[outByte] = (unpacked[i]) |
                                  (unpacked[i + 1] << 3) |
                                  (unpacked[i + 2] << 6);
                packed[outByte + 1] = (unpacked[i + 2] >> 2);
                break;
            case 2:
                packed[outByte] = (unpacked[i]) |
                                  (unpacked[i + 1] << 3);
                break;
            case 1:
                packed[outByte] = (unpacked[i]);
                break;
        }
        return packed;
    }

    template<InternalType INTERNAL_IN,
            ProtoType PROTO_IN,
            InternalType INTERNAL_OUT>
    Sequence<INTERNAL_OUT> packNUMS4(const SequenceProto<INTERNAL_IN, PROTO_IN> &unpacked,
                                     const Alphabet &alphabet) {
        lensq outByte = 0;
        lensq i = 0;
        Sequence<INTERNAL_OUT> packed = reserveSpaceForPacked<INTERNAL_IN, PROTO_IN, INTERNAL_OUT>(unpacked, alphabet);
        for (; i + 8 <= unpacked.size(); i += 8) {
            packed[outByte] = (unpacked[i]) |
                              (unpacked[i + 1] << 4);
            packed[outByte + 1] = (unpacked[i + 2]) |
                                  (unpacked[i + 3] << 4);
            packed[outByte + 2] = (unpacked[i + 4]) |
                                  (unpacked[i + 5] << 4);
            packed[outByte + 3] = (unpacked[i + 6]) |
                                  (unpacked[i + 7] << 4);
            outByte += 4;
        }
        switch (unpacked.size() - i) {
            case 7:
                packed[outByte] = (unpacked[i]) |
                                  (unpacked[i + 1] << 4);
                packed[outByte + 1] = (unpacked[i + 2]) |
                                      (unpacked[i + 3] << 4);
                packed[outByte + 2] = (unpacked[i + 4]) |
                                      (unpacked[i + 5] << 4);
                packed[outByte + 3] = (unpacked[i + 6]);
                break;
            case 6:
                packed[outByte] = (unpacked[i]) |
                                  (unpacked[i + 1] << 4);
                packed[outByte + 1] = (unpacked[i + 2]) |
                                      (unpacked[i + 3] << 4);
                packed[outByte + 2] = (unpacked[i + 4]) |
                                      (unpacked[i + 5] << 4);
                break;
            case 5:
                packed[outByte] = (unpacked[i]) |
                                  (unpacked[i + 1] << 4);
                packed[outByte + 1] = (unpacked[i + 2]) |
                                      (unpacked[i + 3] << 4);
                packed[outByte + 2] = (unpacked[i + 4]);
                break;
            case 4:
                packed[outByte] = (unpacked[i]) |
                                  (unpacked[i + 1] << 4);
                packed[outByte + 1] = (unpacked[i + 2]) |
                                      (unpacked[i + 3] << 4);
                break;
            case 3:
                packed[outByte] = (unpacked[i]) |
                                  (unpacked[i + 1] << 4);
                packed[outByte + 1] = (unpacked[i + 2]);
                break;
            case 2:
                packed[outByte] = (unpacked[i]) |
                                  (unpacked[i + 1] << 4);
                break;
            case 1:
                packed[outByte] = (unpacked[i]);
                break;
        }
        return packed;
    }

    template<InternalType INTERNAL_IN,
            ProtoType PROTO_IN,
            InternalType INTERNAL_OUT>
    Sequence<INTERNAL_OUT> packNUMS5(const SequenceProto<INTERNAL_IN, PROTO_IN> &unpacked,
                                     const Alphabet &alphabet) {
        lensq outByte = 0;
        lensq i = 0;
        Sequence<INTERNAL_OUT> packed = reserveSpaceForPacked<INTERNAL_IN, PROTO_IN, INTERNAL_OUT>(unpacked, alphabet);
        for (; i + 8 <= unpacked.size(); i += 8) {
            packed[outByte] = (unpacked[i]) |
                              (unpacked[i + 1] << 5);
            packed[outByte + 1] = (unpacked[i + 1] >> 3) |
                                  (unpacked[i + 2] << 2) |
                                  (unpacked[i + 3] << 7);
            packed[outByte + 2] = (unpacked[i + 3] >> 1) |
                                  (unpacked[i + 4] << 4);
            packed[outByte + 3] = (unpacked[i + 4] >> 4) |
                                  (unpacked[i + 5] << 1) |
                                  (unpacked[i + 6] << 6);
            packed[outByte + 4] = (unpacked[i + 6] >> 2) |
                                  (unpacked[i + 7] << 3);
            outByte += 5;
        }
        switch (unpacked.size() - i) {
            case 7:
                packed[outByte] = (unpacked[i]) |
                                  (unpacked[i + 1] << 5);
                packed[outByte + 1] = (unpacked[i + 1] >> 3) |
                                      (unpacked[i + 2] << 2) |
                                      (unpacked[i + 3] << 7);
                packed[outByte + 2] = (unpacked[i + 3] >> 1) |
                                      (unpacked[i + 4] << 4);
                packed[outByte + 3] = (unpacked[i + 4] >> 4) |
                                      (unpacked[i + 5] << 1) |
                                      (unpacked[i + 6] << 6);
                packed[outByte + 4] = (unpacked[i + 6] >> 2);
                break;
            case 6:
                packed[outByte] = (unpacked[i]) |
                                  (unpacked[i + 1] << 5);
                packed[outByte + 1] = (unpacked[i + 1] >> 3) |
                                      (unpacked[i + 2] << 2) |
                                      (unpacked[i + 3] << 7);
                packed[outByte + 2] = (unpacked[i + 3] >> 1) |
                                      (unpacked[i + 4] << 4);
                packed[outByte + 3] = (unpacked[i + 4] >> 4) |
                                      (unpacked[i + 5] << 1);
                break;
            case 5:
                packed[outByte] = (unpacked[i]) |
                                  (unpacked[i + 1] << 5);
                packed[outByte + 1] = (unpacked[i + 1] >> 3) |
                                      (unpacked[i + 2] << 2) |
                                      (unpacked[i + 3] << 7);
                packed[outByte + 2] = (unpacked[i + 3] >> 1) |
                                      (unpacked[i + 4] << 4);
                packed[outByte + 3] = (unpacked[i + 4] >> 4);
                break;
            case 4:
                packed[outByte] = (unpacked[i]) |
                                  (unpacked[i + 1] << 5);
                packed[outByte + 1] = (unpacked[i + 1] >> 3) |
                                      (unpacked[i + 2] << 2) |
                                      (unpacked[i + 3] << 7);
                packed[outByte + 2] = (unpacked[i + 3] >> 1);
                break;
            case 3:
                packed[outByte] = (unpacked[i]) |
                                  (unpacked[i + 1] << 5);
                packed[outByte + 1] = (unpacked[i + 1] >> 3) |
                                      (unpacked[i + 2] << 2);
                break;
            case 2:
                packed[outByte] = (unpacked[i]) |
                                  (unpacked[i + 1] << 5);
                packed[outByte + 1] = (unpacked[i + 1] >> 3);
                break;
            case 1:
                packed[outByte] = (unpacked[i]);
                break;
        }
        return packed;
    }


    template<InternalType INTERNAL_IN,
            ProtoType PROTO_IN,
            InternalType INTERNAL_OUT>
    Sequence<INTERNAL_OUT> packNUMS(const SequenceProto<INTERNAL_IN, PROTO_IN> &unpacked,
                                    const Alphabet &alphabet) {
        switch (alphabet.alphabetSize()) {
            case 2:
                return packNUMS2<INTERNAL_IN, PROTO_IN, INTERNAL_OUT>(unpacked, alphabet);
            case 3:
                return packNUMS3<INTERNAL_IN, PROTO_IN, INTERNAL_OUT>(unpacked, alphabet);
            case 4:
                return packNUMS4<INTERNAL_IN, PROTO_IN, INTERNAL_OUT>(unpacked, alphabet);
            case 5:
                return packNUMS5<INTERNAL_IN, PROTO_IN, INTERNAL_OUT>(unpacked, alphabet);
            default:
                throw std::invalid_argument("\"alphabet\" has bad alphabet size");
        }
    }
}

#endif //TIDYSQ_PACKNUMS_H
