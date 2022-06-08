#ifndef POSCLASS_H
#define POSCLASS_H

#include <string>
#include <experimental/string_view>
#include <map>
#include <vector>
#include <iostream>

template<class T>
struct PosclassImpl
{
    int chr;

    int position;

    std::string ref;

    std::string alt;

    PosclassImpl(): chr(0), position(-1), ref(""), alt(""){}

    static const std::map<std::string, std::string>& Convmap() {static const std::map<std::string, std::string> convmap={{"A","T"},{"C","G"}, {"G","C"}, {"T","A"}, {"-","-"}, {"I","I"}, {"D","D"}};return convmap;}

    std::string CreateStringMap(const std::string& ref, const std::string& alt) const
    {
        return static_cast<const T*>(this)->CreateStringMapImpl(ref, alt);
    }

    std::tuple<int,int,std::string> CreateTuplePosition() const
    {
        auto refalt = CreateStringMap(ref, alt);

        int pos = position;

        if(ref.size() > 1 || alt.size()>1)
        {
            pos += 1;

            refalt = "0," + refalt;
        }

        return std::make_tuple(chr, position, refalt);
    }

    PosclassImpl(int chr_, int position_, const std::string& ref_, const std::string& alt_): chr(chr_), position(position_), ref(ref_), alt(alt_)
    {
    }

    PosclassImpl(int chr_, int position_, const std::experimental::string_view& ref_, const std::experimental::string_view& alt_): chr(chr_), position(position_), ref(ref_), alt(alt_)
    {
    }

    PosclassImpl(const SNPclass& snpclass):chr(snpclass.chr), position(snpclass.position), ref(snpclass.ref), alt(snpclass.alt)
    {
    }

    PosclassImpl(const Imputeclass& imputeclass):chr(ConvertChrStr(imputeclass.chromosome)), position(imputeclass.position), ref(imputeclass.ref), alt(imputeclass.alt)
    {
    }

    std::pair<int,int> getPos() const {return {chr,position};}

    bool operator<(const PosclassImpl& rhs) const
    {
        auto pos1 = std::make_pair(chr, position);

        auto pos2 = std::make_pair(rhs.chr, rhs.position);

        if(pos1 != pos2)
        {
            return pos1 < pos2;
        }
        else
        {
            return (CreateStringMap(ref, alt) < CreateStringMap(rhs.ref, rhs.alt));
        }
    }

    bool operator<=(const PosclassImpl& rhs) const
    {
        auto pos1 = std::make_pair(chr, position);

        auto pos2 = std::make_pair(rhs.chr, rhs.position);

        if(pos1 != pos2)
        {
            return pos1 < pos2;
        }
        else
        {
            return (CreateStringMap(ref, alt) <= CreateStringMap(rhs.ref, rhs.alt));
        }
    }

    bool operator>(const PosclassImpl& rhs) const
    {
        auto pos1 = std::make_pair(chr, position);

        auto pos2 = std::make_pair(rhs.chr, rhs.position);

        if(pos1 != pos2)
        {
            return pos1 > pos2;
        }
        else
        {
            return (CreateStringMap(ref, alt) > CreateStringMap(rhs.ref, rhs.alt));
        }
    }

    bool operator==(const PosclassImpl& rhs) const
    {
        auto pos1 = std::make_pair(chr, position);

        auto pos2 = std::make_pair(rhs.chr, rhs.position);

        if(pos1 != pos2)
        {
            return false;
        }
        else
        {
            return (CreateStringMap(ref, alt) == CreateStringMap(rhs.ref, rhs.alt));
        }
    }
};


struct Posclass: PosclassImpl<Posclass>
{
    Posclass(int chr_=-1, int position_=-1, const std::string& ref_="", const std::string& alt_=""): PosclassImpl<Posclass>(chr_,position_,ref_,alt_)
    {
    }

    Posclass(int chr_, int position_, const std::experimental::string_view& ref_, const std::experimental::string_view& alt_): PosclassImpl<Posclass>(chr_,position_,ref_,alt_)
    {
    }

    Posclass(const SNPclass& snpclass):PosclassImpl<Posclass>(snpclass)
    {
    }

    Posclass(const Imputeclass& imputeclass):PosclassImpl<Posclass>(imputeclass)
    {
    }

    std::string CreateStringMapImpl(const std::string& ref, const std::string& alt) const
    {
        std::vector<std::string> refalt(2);

        refalt[0] = ref;
        refalt[1] = alt;

        std::vector<std::string> refalt2(2);

        if(refalt[0].size() == 1 && refalt[1].size() == 1)
        {
            try
            {
                refalt2[0] = Convmap().at(refalt[0]);

                refalt2[1] = Convmap().at(refalt[1]);
            }
            catch(...)
            {
                std::cout << "Mismatch for " << refalt[0] << '\t' << refalt[1] << '\n';
            }

            std::sort(refalt.begin(), refalt.end());

            std::sort(refalt2.begin(), refalt2.end());

            if(refalt > refalt2) refalt = refalt2;

            refalt[0] = "z" + refalt[0];
        }
        else
        {
    //		refalt[0] = (refalt[0].size() > 1 ? "I" : "D");

    //		refalt[1] = (refalt[1].size() > 1 ? "I" : "D");

//            std::sort(refalt.begin(), refalt.end());

// don't sort for insertion deletions as order matters
        }

        return refalt[0] + ',' + refalt[1];
    }
};


struct PosclassID: PosclassImpl<PosclassID>
{
    PosclassID(int chr_=-1, int position_=-1, const std::string& ref_="", const std::string& alt_=""): PosclassImpl<PosclassID>(chr_,position_,ref_,alt_)
    {
    }

    PosclassID(const SNPclass& snpclass):PosclassImpl<PosclassID>(snpclass)
    {
    }

    std::string CreateStringMapImpl(const std::string& ref, const std::string& alt) const
    {
        std::vector<std::string> refalt(2);

        refalt[0] = ref;
        refalt[1] = alt;

        std::vector<std::string> refalt2(2);

        if(refalt[0].size() == 1 && refalt[1].size() == 1)
        {
            try
            {
                refalt2[0] = Convmap().at(refalt[0]);

                refalt2[1] = Convmap().at(refalt[1]);
            }
            catch(...)
            {
                std::cout << "Mismatch for " << refalt[0] << '\t' << refalt[1] << '\n';
            }

            std::sort(refalt.begin(), refalt.end());

            std::sort(refalt2.begin(), refalt2.end());

            if(refalt > refalt2) refalt = refalt2;

            refalt[0] = "z" + refalt[0];
        }
        else
        {
            refalt[0] = (refalt[0].size() > refalt[1].size() ? "I" : "D");

            refalt[1] = (refalt[1].size() > refalt[0].size() ? "I" : "D");

//            std::sort(refalt.begin(), refalt.end());
        }

        return refalt[0] + ',' + refalt[1];
    }
};


struct PosclassUnordered: PosclassImpl<PosclassUnordered>
{
    PosclassUnordered(int chr_=-1, int position_=-1, const std::string& ref_="", const std::string& alt_=""): PosclassImpl<PosclassUnordered>(chr_,position_,ref_,alt_)
    {
    }

    PosclassUnordered(const SNPclass& snpclass):PosclassImpl<PosclassUnordered>(snpclass)
    {
    }

    std::string CreateStringMapImpl(const std::string& ref, const std::string& alt) const
    {
        std::vector<std::string> refalt(2);

        refalt[0] = ref;
        refalt[1] = alt;

        std::vector<std::string> refalt2(2);

        if(refalt[0].size() == 1 && refalt[1].size() == 1)
        {
            try
            {
                refalt2[0] = Convmap().at(refalt[0]);

                refalt2[1] = Convmap().at(refalt[1]);
            }
            catch(...)
            {
//                std::cout << "Mismatch for " << refalt[0] << '\t' << refalt[1] << '\n';
            }

            std::sort(refalt.begin(), refalt.end());

            std::sort(refalt2.begin(), refalt2.end());

            if(refalt > refalt2) refalt = refalt2;

            refalt[0] = "z" + refalt[0];
        }
        else
        {
            std::sort(refalt.begin(), refalt.end());
        }

        return refalt[0] + ',' + refalt[1];
    }
};


namespace
{

template<class T>
std::istream& operator>>(std::istream& is, PosclassImpl<T>& v)
{
    is >> v.chr >> v.position >> v.ref >> v.alt;

    return is;
}

template<class T>
std::ostream& operator<<(std::ostream& os, const PosclassImpl<T>& v)
{
    os << v.chr << '\t' << v.position << '\t' << v.ref << '\t' << v.alt;

    return os;
}

}

#endif // POSCLASS_H



