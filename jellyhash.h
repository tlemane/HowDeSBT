// jellyhash.h-- an implementation of jellyfish's hash function, using a hard-
//   coded matrix and supporting only 20-mers.
//   (based on code in jellyfish's mer_dna.hpp and rectangular_binary_matrix.hpp)
//
// Jellyfish's hash function is a matrix multiply in GF(2).  The matrix herein
// is *a* matrix generated by jellyfish for 20-mers.

#ifndef jellyhash_H
#define jellyhash_H

#define jellyhash_K     20

#define jellyhash_Abits 0  // A in jellyfish/mer_dna.hpp
#define jellyhash_Cbits 1  // C in jellyfish/mer_dna.hpp
#define jellyhash_Gbits 2  // G in jellyfish/mer_dna.hpp
#define jellyhash_Tbits 3  // T in jellyfish/mer_dna.hpp

#define jellyhash_NonZero 2*1000*1000*1000


static const std::uint64_t jellyMatrix[2*jellyhash_K] =
	{
	0x9C9EC8F1EB8B4567, 0xED2C371426334873, 0x8895563B2AE8944A, 0xBF46D62EC6E87CCD,
	0x206DC7BEEEB141F2, 0x3457C01F7545E146, 0xA36C49FE12200854, 0xD4643379DF16E9E8,
	0xECD48956940E0F76, 0xEFF7308CCDED7263, 0x7ADA01E6C1A7C4C9, 0x6466D40365E45D32,
	0x2F20F9163F2DBA31, 0xC0DB184922BBD95A, 0x9C8768C5F33AB105, 0x59D8E117AD1D5AE9,
	0xCE60E32C88EDBDAB, 0x46269DA6CB03E0C6, 0xEB2A218431F32454, 0x1EA57E1D02901D82,
	0xDF0F6F4F5E7FF521, 0x9888869C6CEAF087, 0x0853F5287006C83E, 0x6102EB7F5577F8E1,
	0xFDD197C07804823E, 0xC918EE7ADC482A97, 0xFB545DE5D1EAD36B, 0xAE155617153EA438,
	0xBA921F2C6A2342EC, 0xFB36268CF25A06FB, 0xD2D623D57A6D8D3C, 0xCE10DFF6ADE91B18,
	0xFA129206B2FFF902, 0xAF70903EB49ABB43, 0xAD71B0CEF9A1DEAA, 0x7483BB7430C6A529,
	0xC8FE704F4F4EF005, 0x8E4E195D675AC794, 0x98D6EA18580115BE, 0xB56D6BD7354FE9F9
	};


class JellyHash
	{
public:
	unsigned int   k;
	unsigned int   numMerBits;
	unsigned int   numMerWords;
	std::uint64_t  nucToBits[256];
	std::uint64_t  hForward;

public:
	JellyHash
	   (const unsigned int _k,
		const std::uint64_t _seed=0,
		const bool _allowN=false)
	  :	k(_k)
		{
		if (_k != jellyhash_K)
			fatal ("internal error: JellyHash(" + std::to_string(_k) + ") is not supported");

		numMerBits  = 2*_k;
		numMerWords = (numMerBits + 63) / 64;

		for (unsigned int ch=0 ; ch<256 ; ch++)
			nucToBits[ch] = (std::uint64_t) -1;

		nucToBits    [(unsigned char)'A']
		  = nucToBits[(unsigned char)'a'] = jellyhash_Abits;

		nucToBits    [(unsigned char)'C']
		  = nucToBits[(unsigned char)'c'] = jellyhash_Cbits;

		nucToBits    [(unsigned char)'G'] 
		  = nucToBits[(unsigned char)'g'] = jellyhash_Gbits;

		nucToBits    [(unsigned char)'T'] 
		  = nucToBits[(unsigned char)'t'] = jellyhash_Tbits;
		}

	inline std::uint64_t hash
	   (const char* s)  // only first k characters are used
		{
		// note that we return 0 iff we don't have k valid characters
		std::uint64_t data[numMerWords];

		for (unsigned int merIx=0 ; merIx<numMerWords ; merIx++)
			data[merIx] = 0;

		unsigned int bitPos, merIx;
		for (unsigned int ix=0 ; ix<k ; ix++)
			{
			std::uint64_t twoBits = nucToBits[(unsigned char) s[ix]];
			if (twoBits == (std::uint64_t) -1) // chIn not in {A,C,G,T}
				return 0;

			bitPos = 2*((k-1)-ix);
			merIx  = bitPos / 64;
			bitPos = bitPos - (64*merIx);
			data[merIx] |= twoBits << bitPos;
			}

		hForward = apply_matrix(data,k);
		return (hForward == 0)? jellyhash_NonZero : hForward;  // differs from jellyfish
		}

	inline std::uint64_t hash
	   (const std::string& s)  // only first k characters are used
		{
		// note that we return 0 iff we don't have k valid characters
		std::uint64_t data[numMerWords];

		if (s.length() < k) return 0;

		for (unsigned int merIx=0 ; merIx<numMerWords ; merIx++)
			data[merIx] = 0;

		unsigned int bitPos, merIx;
		for (unsigned int ix=0 ; ix<k ; ix++)
			{
			std::uint64_t twoBits = nucToBits[(unsigned char) s[ix]];
			if (twoBits == (std::uint64_t) -1) // chIn not in {A,C,G,T}
				return 0;

			bitPos = 2*((k-1)-ix);
			merIx  = bitPos / 64;
			bitPos = bitPos - (64*merIx);
			data[merIx] |= twoBits << bitPos;
			}

		hForward = apply_matrix(data,k);
		return (hForward == 0)? jellyhash_NonZero : hForward;  // differs from jellyfish
		}

	inline std::uint64_t hash
	   (const std::uint64_t* _data)
		{
		hForward = apply_matrix(_data,k);
		return (hForward == 0)? jellyhash_NonZero : hForward;  // differs from jellyfish
		}

	static inline void fill_hash_values
	   (std::uint64_t       hashValues[],
		const int           numHashes,
		const std::uint64_t h1,
		const std::uint64_t h2)
		{
		fatal ("internal error: JellyHash::fill_hash_values() is not implemented");
		}

	// see RectangularBinaryMatrix::times_loop in jellyfish/rectangular_binary_matrix.hpp
	static inline std::uint64_t apply_matrix
	   (const std::uint64_t* _data,
		const unsigned int   _k)
		{
		unsigned int   numMerBits  = 2*_k;
		unsigned int   numMerWords = (numMerBits + 63) / 64;
		unsigned int   extraBits   = numMerBits & 0x3F;
		std::uint64_t* data = (std::uint64_t*) _data;
		std::uint64_t* m    = (std::uint64_t*) &jellyMatrix[numMerBits-1];
		std::uint64_t  h    = 0;
		std::uint64_t  x    = 0;

		std::uint64_t bitsInWord = 64;
		for (unsigned int i=0 ; i<numMerWords ; i++)
			{
			x = data[i];
			if ((i < numMerWords-1) || (extraBits == 0))
				bitsInWord = 64;
			else
				{
				bitsInWord = extraBits;
				x &= (((std::uint64_t)1)<<bitsInWord)-1;    // (not really necessary)
				}

			for ( ; bitsInWord>=8 ; bitsInWord-=8,m-=8)
				{
				// nota bene:
				//    lsbit(x)  x&1   -(x&1)
				//       0       0     0, all zeros
				//       1       1    -1, all ones
				h ^= (-(x&1)) & m[ 0]; x >>= 1;
				h ^= (-(x&1)) & m[-1]; x >>= 1;
				h ^= (-(x&1)) & m[-2]; x >>= 1;
				h ^= (-(x&1)) & m[-3]; x >>= 1;
				h ^= (-(x&1)) & m[-4]; x >>= 1;
				h ^= (-(x&1)) & m[-5]; x >>= 1;
				h ^= (-(x&1)) & m[-6]; x >>= 1;
				h ^= (-(x&1)) & m[-7]; x >>= 1;
				}
			}

		switch (bitsInWord)
			{
			case 7: h ^= (-(x&1)) & *(m--); x >>= 1;
			case 6: h ^= (-(x&1)) & *(m--); x >>= 1;
			case 5: h ^= (-(x&1)) & *(m--); x >>= 1;
			case 4: h ^= (-(x&1)) & *(m--); x >>= 1;
			case 3: h ^= (-(x&1)) & *(m--); x >>= 1;
			case 2: h ^= (-(x&1)) & *(m--); x >>= 1;
			case 1: h ^= (-(x&1)) & *(m--);
			}

		return h;
		}
	};

class JellyHashCanonical
	{
public:
	unsigned int   k;
	unsigned int   numMerBits;
	unsigned int   numMerWords;
	std::uint64_t  nucToBits[256], nucToBitsRC[256], revCompBits[4];
	std::uint64_t  hCanonical;

public:
	JellyHashCanonical
	   (const unsigned int _k,
		const std::uint64_t _seed=0,
		const bool _allowN=false)
	  :	k(_k)
		{
		if (_k != jellyhash_K)
			fatal ("internal error: JellyHashCanonical(" + std::to_string(_k) + ") is not supported");

		numMerBits  = 2*_k;
		numMerWords = (numMerBits + 63) / 64;

		for (unsigned int ch=0 ; ch<256 ; ch++)
			nucToBits[ch] = (std::uint64_t) -1;

		nucToBits    [(unsigned char)'A']
		  = nucToBits[(unsigned char)'a'] = jellyhash_Abits;

		nucToBits    [(unsigned char)'C']
		  = nucToBits[(unsigned char)'c'] = jellyhash_Cbits;

		nucToBits    [(unsigned char)'G'] 
		  = nucToBits[(unsigned char)'g'] = jellyhash_Gbits;

		nucToBits    [(unsigned char)'T'] 
		  = nucToBits[(unsigned char)'t'] = jellyhash_Tbits;

		nucToBitsRC    [(unsigned char)'A'] 
		  = nucToBitsRC[(unsigned char)'a'] = jellyhash_Tbits;

		nucToBitsRC    [(unsigned char)'C'] 
		  = nucToBitsRC[(unsigned char)'c'] = jellyhash_Gbits;

		nucToBitsRC    [(unsigned char)'G'] 
		  = nucToBitsRC[(unsigned char)'g'] = jellyhash_Cbits;

		nucToBitsRC    [(unsigned char)'T']
		  = nucToBitsRC[(unsigned char)'t'] = jellyhash_Abits;

		revCompBits[jellyhash_Abits] = jellyhash_Tbits;
		revCompBits[jellyhash_Cbits] = jellyhash_Gbits;
		revCompBits[jellyhash_Gbits] = jellyhash_Cbits;
		revCompBits[jellyhash_Tbits] = jellyhash_Abits;
		}

	inline std::uint64_t hash
	   (const char* s)  // only first k characters are used
		{
		// note that we return 0 iff we don't have k valid characters
		std::uint64_t data[numMerWords], dataRC[numMerWords];

		for (unsigned int merIx=0 ; merIx<numMerWords ; merIx++)
			data[merIx] = dataRC[merIx] = 0;

		unsigned int bitPos, merIx, bitPosRC, merIxRC;
		for (unsigned int ix=0 ; ix<k ; ix++)
			{
			std::uint64_t twoBits = nucToBits[(unsigned char) s[ix]];
			if (twoBits == (std::uint64_t) -1) // chIn not in {A,C,G,T}
				return 0;

			bitPos   = 2*((k-1)-ix);
			merIx    = bitPos / 64;
			bitPos   = bitPos - (64*merIx);

			bitPosRC = 2*ix;
			merIxRC  = bitPosRC / 64;
			bitPosRC = bitPosRC - (64*merIxRC);

			data  [merIx  ] |= twoBits << bitPos;
			dataRC[merIxRC] |= revCompBits[twoBits] << bitPosRC;
			}

		if (dataRC[0] < data[0])  // $$$ only correct for k <= 32
			{
			for (unsigned int merIx=0 ; merIx<numMerWords ; merIx++)
				data[merIx] = dataRC[merIx];
			}

		hCanonical = JellyHash::apply_matrix(data,k);
		return (hCanonical == 0)? jellyhash_NonZero : hCanonical;  // differs from jellyfish
		}

	inline std::uint64_t hash
	   (const std::string& s)  // only first k characters are used
		{
		// note that we return 0 iff we don't have k valid characters
		std::uint64_t data[numMerWords], dataRC[numMerWords];

		if (s.length() < k) return 0;

		for (unsigned int merIx=0 ; merIx<numMerWords ; merIx++)
			data[merIx] = dataRC[merIx] = 0;

		unsigned int bitPos, merIx, bitPosRC, merIxRC;
		for (unsigned int ix=0 ; ix<k ; ix++)
			{
			std::uint64_t twoBits = nucToBits[(unsigned char) s[ix]];
			if (twoBits == (std::uint64_t) -1) // chIn not in {A,C,G,T}
				return 0;

			bitPos   = 2*((k-1)-ix);
			merIx    = bitPos / 64;
			bitPos   = bitPos - (64*merIx);

			bitPosRC = 2*ix;
			merIxRC  = bitPosRC / 64;
			bitPosRC = bitPosRC - (64*merIxRC);

			data  [merIx  ] |= twoBits << bitPos;
			dataRC[merIxRC] |= revCompBits[twoBits] << bitPosRC;
			}


		if (dataRC[0] < data[0])  // $$$ only correct for k <= 32
			{
			for (unsigned int merIx=0 ; merIx<numMerWords ; merIx++)
				data[merIx] = dataRC[merIx];
			}

		hCanonical = JellyHash::apply_matrix(data,k);
		return (hCanonical == 0)? jellyhash_NonZero : hCanonical;  // differs from jellyfish
		}

	inline std::uint64_t hash
	   (const std::uint64_t* _data)
		{
		std::uint64_t dataRC[numMerWords], dataCanon[numMerWords];

		for (unsigned int merIx=0 ; merIx<numMerWords ; merIx++)
			dataRC[merIx] = 0;

		unsigned int bitPos, merIx, bitPosRC, merIxRC;
		for (unsigned int ix=0 ; ix<k ; ix++)
			{
			bitPos   = 2*((k-1)-ix);
			merIx    = bitPos / 64;
			bitPos   = bitPos - (64*merIx);

			bitPosRC = 2*ix;
			merIxRC  = bitPosRC / 64;
			bitPosRC = bitPosRC - (64*merIxRC);

			std::uint64_t twoBits = (_data[merIx]>>bitPos) & 3;
			dataRC[merIxRC] |= revCompBits[twoBits] << bitPosRC;
			}

		if (dataRC[0] < _data[0])  // $$$ only correct for k <= 32
			{
			for (unsigned int merIx=0 ; merIx<numMerWords ; merIx++)
				dataCanon[merIx] = dataRC[merIx];
			}
		else
			{
			for (unsigned int merIx=0 ; merIx<numMerWords ; merIx++)
				dataCanon[merIx] = _data[merIx];
			}

		hCanonical = JellyHash::apply_matrix(dataCanon,k);
		return (hCanonical == 0)? jellyhash_NonZero : hCanonical;  // differs from jellyfish
		}
	};


#undef jellyhash_K
#undef jellyhash_Abits
#undef jellyhash_Cbits
#undef jellyhash_Gbits
#undef jellyhash_Tbits
#undef jellyhash_NonZero


#endif // jellyhash_H