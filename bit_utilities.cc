// bit_utilities.cc-- bit-related utility functions.

#include <cstdint>

#include "bit_utilities.h"

#define u8  std::uint8_t
#define u64 std::uint64_t

#define least_significant(type,numBits) ((((type)1)<<(numBits))-1)

//----------
//
// lookup table(s)
//
//----------

// popCount8 maps a byte to the number of 1s in that byte

static const u8 popCount8[256] =
	{
	0,1,1,2,1,2,2,3,1,2,2,3,2,3,3,4,
	1,2,2,3,2,3,3,4,2,3,3,4,3,4,4,5,
	1,2,2,3,2,3,3,4,2,3,3,4,3,4,4,5,
	2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,
	1,2,2,3,2,3,3,4,2,3,3,4,3,4,4,5,
	2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,
	2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,
	3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,
	1,2,2,3,2,3,3,4,2,3,3,4,3,4,4,5,
	2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,
	2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,
	3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,
	2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,
	3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,
	3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,
	4,5,5,6,5,6,6,7,5,6,6,7,6,7,7,8
	};

//----------
//
// bitwise_and, bitwise_mask, bitwise_or, bitwise_or_not, bitwise_xor,
// and bitwise_xnor--
//	Create the bitwise AND, ANDNOT, OR, XOR, or XNOR of two bit arrays.
//
//----------
//
// Arguments:
//	const void*	bits1, bits2:	Bit arrays to read.
//	void*		dstBits:		Bit array to fill.
//	u64			numBits:		The length of the bit arrays, counted in *bits*.
//								.. See note (1) below.
//
// Returns:
//	(nothing)
//
//----------
//
// Notes:
//	(1)	The number of bytes in the bit arrays is ceil(numBits/8). When numBits
//		is not a multiple of 8, the remaining bits are read from the least
//		significant bits of the final byte.
//	(2)	We process the bytes in 64-bit chunks until we get to the final chunk.
//		The final chunk is processed byte-by-byte, so that we do not access
//		any bytes beyond the bit arrays.
//	(3)	Equivalences to set operations:
//		  function        logic operation   set operation
//		  ------------    ---------------   --------------------
//		  bitwise_and     a AND b           a intersect b
//		  bitwise_mask    a AND (NOT b)     a\b (set difference)
//		  bitwise_or      a OR b            a union b
//		  bitwise_or_not  a OR (NOT b)      a union complement(b)
//		  bitwise_xor     a XOR b           a\b union b\a
//		  bitwise_xnor    NOT (a XOR b)     a == b
//
//----------

void bitwise_and
   (const void*	bits1,
	const void*	bits2,
	void*		dstBits,
	const u64	numBits)
	{
	u64*		scan1 = (u64*) bits1;
	u64*		scan2 = (u64*) bits2;
	u64*		dst   = (u64*) dstBits;
	u64			n;

	for (n=numBits ; n>=64 ; n-=64)
		*(dst++) = *(scan1++) & *(scan2++);

	if (n == 0) return;

	u8*	scan1b = (u8*) scan1;
	u8*	scan2b = (u8*) scan2;
	u8*	dstb   = (u8*) dst;
	for ( ; n>=8 ; n-=8)
		*(dstb++) = *(scan1b++) & *(scan2b++);

	if (n == 0) return;

	u8 mask = least_significant(u8,n);
	*dstb = (*scan1b & *scan2b) & mask;  // leftover bits intentionally set to zero
	}


void bitwise_and
   (void*		dstBits,
	const void*	bits2,
	const u64	numBits)
	{
	u64*		dst   = (u64*) dstBits;
	u64*		scan2 = (u64*) bits2;
	u64			n;

	for (n=numBits ; n>=64 ; n-=64,dst++)
		*dst &= *(scan2++);

	if (n == 0) return;

	u8*	scan2b = (u8*) scan2;
	u8*	dstb   = (u8*) dst;
	for ( ; n>=8 ; n-=8,dstb++)
		*dstb &= *(scan2b++);

	if (n == 0) return;

	u8 mask = least_significant(u8,n);
	*dstb = ((*dstb & *scan2b) & mask) | (*dstb & ~mask);
	}


void bitwise_mask
   (const void*	bits1,
	const void*	bits2,
	void*		dstBits,
	const u64	numBits)
	{
	u64*		scan1 = (u64*) bits1;
	u64*		scan2 = (u64*) bits2;
	u64*		dst   = (u64*) dstBits;
	u64			n;

	for (n=numBits ; n>=64 ; n-=64)
		*(dst++) = *(scan1++) & ~*(scan2++);

	if (n == 0) return;

	u8*	scan1b = (u8*) scan1;
	u8*	scan2b = (u8*) scan2;
	u8*	dstb   = (u8*) dst;
	for ( ; n>=8 ; n-=8)
		*(dstb++) = *(scan1b++) & ~*(scan2b++);

	if (n == 0) return;

	u8 mask = least_significant(u8,n);
	*dstb = (*scan1b & ~*scan2b) & mask;  // leftover bits intentionally set to zero
	}


void bitwise_mask
   (void*		dstBits,
	const void*	bits2,
	const u64	numBits)
	{
	u64*		dst   = (u64*) dstBits;
	u64*		scan2 = (u64*) bits2;
	u64			n;

	for (n=numBits ; n>=64 ; n-=64,dst++)
		*dst &= ~*(scan2++);

	if (n == 0) return;

	u8*	scan2b = (u8*) scan2;
	u8*	dstb   = (u8*) dst;
	for ( ; n>=8 ; n-=8,dstb++)
		*dstb &= ~*(scan2b++);

	if (n == 0) return;

	u8 mask = least_significant(u8,n);
	*dstb = ((*dstb & ~*scan2b) & mask) | (*dstb & ~mask);
	}


void bitwise_or
   (const void*	bits1,
	const void*	bits2,
	void*		dstBits,
	const u64	numBits)
	{
	u64*		scan1 = (u64*) bits1;
	u64*		scan2 = (u64*) bits2;
	u64*		dst   = (u64*) dstBits;
	u64			n;

	for (n=numBits ; n>=64 ; n-=64)
		*(dst++) = *(scan1++) | *(scan2++);

	if (n == 0) return;

	u8*	scan1b = (u8*) scan1;
	u8*	scan2b = (u8*) scan2;
	u8*	dstb   = (u8*) dst;
	for ( ; n>=8 ; n-=8)
		*(dstb++) = *(scan1b++) | *(scan2b++);

	if (n == 0) return;

	u8 mask = least_significant(u8,n);
	*dstb = (*scan1b | *scan2b) & mask;  // leftover bits intentionally set to zero
	}


void bitwise_or
   (void*		dstBits,
	const void*	bits2,
	const u64	numBits)
	{
	u64*		dst   = (u64*) dstBits;
	u64*		scan2 = (u64*) bits2;
	u64			n;

	for (n=numBits ; n>=64 ; n-=64,dst++)
		*dst |= *(scan2++);

	if (n == 0) return;

	u8*	scan2b = (u8*) scan2;
	u8*	dstb   = (u8*) dst;
	for ( ; n>=8 ; n-=8,dstb++)
		*dstb |= *(scan2b++);

	if (n == 0) return;

	u8 mask = least_significant(u8,n);
	*dstb = ((*dstb | *scan2b) & mask) | (*dstb & ~mask);
	}


void bitwise_or_not
   (const void*	bits1,
	const void*	bits2,
	void*		dstBits,
	const u64	numBits)
	{
	u64*		scan1 = (u64*) bits1;
	u64*		scan2 = (u64*) bits2;
	u64*		dst   = (u64*) dstBits;
	u64			n;

	for (n=numBits ; n>=64 ; n-=64)
		*(dst++) = *(scan1++) | ~*(scan2++);

	if (n == 0) return;

	u8*	scan1b = (u8*) scan1;
	u8*	scan2b = (u8*) scan2;
	u8*	dstb   = (u8*) dst;
	for ( ; n>=8 ; n-=8)
		*(dstb++) = *(scan1b++) | ~*(scan2b++);

	if (n == 0) return;

	u8 mask = least_significant(u8,n);
	*dstb = (*scan1b | ~*scan2b) & mask;  // leftover bits intentionally set to zero
	}


void bitwise_or_not
   (void*		dstBits,
	const void*	bits2,
	const u64	numBits)
	{
	u64*		dst   = (u64*) dstBits;
	u64*		scan2 = (u64*) bits2;
	u64			n;

	for (n=numBits ; n>=64 ; n-=64,dst++)
		*dst |= ~*(scan2++);

	if (n == 0) return;

	u8*	scan2b = (u8*) scan2;
	u8*	dstb   = (u8*) dst;
	for ( ; n>=8 ; n-=8,dstb++)
		*dstb |= ~*(scan2b++);

	if (n == 0) return;

	u8 mask = least_significant(u8,n);
	*dstb = ((*dstb | ~*scan2b) & mask) | (*dstb & ~mask);
	}


void bitwise_xor
   (const void*	bits1,
	const void*	bits2,
	void*		dstBits,
	const u64	numBits)
	{
	u64*		scan1 = (u64*) bits1;
	u64*		scan2 = (u64*) bits2;
	u64*		dst   = (u64*) dstBits;
	u64			n;

	for (n=numBits ; n>=64 ; n-=64)
		*(dst++) = *(scan1++) ^ *(scan2++);

	if (n == 0) return;

	u8*	scan1b = (u8*) scan1;
	u8*	scan2b = (u8*) scan2;
	u8*	dstb   = (u8*) dst;
	for ( ; n>=8 ; n-=8)
		*(dstb++) = *(scan1b++) ^ *(scan2b++);

	if (n == 0) return;

	u8 mask = least_significant(u8,n);
	*dstb = (*scan1b ^ *scan2b) & mask;  // leftover bits intentionally set to zero
	}


void bitwise_xor
   (void*		dstBits,
	const void*	bits2,
	const u64	numBits)
	{
	u64*		dst   = (u64*) dstBits;
	u64*		scan2 = (u64*) bits2;
	u64			n;

	for (n=numBits ; n>=64 ; n-=64,dst++)
		*dst ^= *(scan2++);

	if (n == 0) return;

	u8*	scan2b = (u8*) scan2;
	u8*	dstb   = (u8*) dst;
	for ( ; n>=8 ; n-=8,dstb++)
		*dstb ^= *(scan2b++);

	if (n == 0) return;

	u8 mask = least_significant(u8,n);
	*dstb = ((*dstb ^ *scan2b) & mask) | (*dstb & ~mask);
	}


void bitwise_xnor
   (const void*	bits1,
	const void*	bits2,
	void*		dstBits,
	const u64	numBits)
	{
	u64*		scan1 = (u64*) bits1;
	u64*		scan2 = (u64*) bits2;
	u64*		dst   = (u64*) dstBits;
	u64			n;

	for (n=numBits ; n>=64 ; n-=64)
		*(dst++) = ~(*(scan1++) ^ *(scan2++));

	if (n == 0) return;

	u8*	scan1b = (u8*) scan1;
	u8*	scan2b = (u8*) scan2;
	u8*	dstb   = (u8*) dst;
	for ( ; n>=8 ; n-=8)
		*(dstb++) = ~(*(scan1b++) ^ *(scan2b++));

	if (n == 0) return;

	u8 mask = least_significant(u8,n);
	*dstb = ~((*scan1b ^ *scan2b)) & mask;  // leftover bits intentionally set to zero
	}


void bitwise_xnor
   (void*		dstBits,
	const void*	bits2,
	const u64	numBits)
	{
	u64*		dst   = (u64*) dstBits;
	u64*		scan2 = (u64*) bits2;
	u64			n;

	for (n=numBits ; n>=64 ; n-=64,dst++)
		*dst = ~(*dst ^ *(scan2++));

	if (n == 0) return;

	u8*	scan2b = (u8*) scan2;
	u8*	dstb   = (u8*) dst;
	for ( ; n>=8 ; n-=8,dstb++)
		*dstb = ~(*dstb ^ *(scan2b++));

	if (n == 0) return;

	u8 mask = least_significant(u8,n);
	*dstb = (~((*dstb ^ *scan2b)) & mask) | (*dstb & ~mask);
	}

//----------
//
// bitwise_complement--
//	Create the bitwise complement of a bit array.
//
//----------
//
// Arguments:
//	const void*	bits:		Bit array to read.
//	void*		dstBits:	Bit array to fill.
//	u64			numBits:	The length of the bit arrays, counted in *bits*.
//							.. See note (1) below.
//
// Returns:
//	(nothing)
//
//----------
//
// Notes:
//	(1)	The number of bytes in the bit arrays is ceil(numBits/8). When numBits
//		is not a multiple of 8, the remaining bits are read from the least
//		significant bits of the final byte.
//	(2)	We process the bytes in 64-bit chunks until we get to the final chunk.
//		The final chunk is processed byte-by-byte, so that we do not access
//		any bytes beyond the bit arrays.
//
//----------

void bitwise_complement
   (const void*	bits,
	void*		dstBits,
	const u64	numBits)
	{
	u64*		scan = (u64*) bits;
	u64*		dst  = (u64*) dstBits;
	u64			n;

	for (n=numBits ; n>=64 ; n-=64)
		*(dst++) = ~*(scan++);

	if (n == 0) return;

	u8*	scanb = (u8*) scan;
	u8*	dstb  = (u8*) dst;
	for ( ; n>=8 ; n-=8)
		*(dstb++) = ~*(scanb++);

	if (n == 0) return;

	u8 mask = least_significant(u8,n);
	*dstb = ~((*scanb) & mask);  // leftover bits intentionally set to zero
	}


void bitwise_complement
   (void*		dstBits,
	const u64	numBits)
	{
	u64*		dst = (u64*) dstBits;
	u64			n;

	for (n=numBits ; n>=64 ; n-=64,dst++)
		*dst = ~*dst;

	if (n == 0) return;

	u8*	dstb = (u8*) dst;
	for ( ; n>=8 ; n-=8,dstb++)
		*dstb = ~*dstb;

	if (n == 0) return;

	u8 mask = least_significant(u8,n);
	*dstb = ~((*dstb) & mask) | (*dstb & ~mask);
	}

//----------
//
// bitwise_fill--
//	Fill a bit array with 1s or 0s.
//
//----------
//
// Arguments:
//	void*		dstBits:	Bit array to fill.
//	int			bitVal:		Value to fill the bit array with.  Only the least
//							.. significant bit is used (the other bits are
//							.. ignored).
//	u64			numBits:	The length of the bit arrays, counted in *bits*.
//							.. See note (1) below.
//
// Returns:
//	(nothing)
//
//----------
//
// Notes:
//	(1)	The number of bytes in the bit arrays is ceil(numBits/8). When numBits
//		is not a multiple of 8, the remaining bits are in the least significant
//		bits of the final byte.
//	(2)	We process the bytes in 64-bit chunks until we get to the final chunk.
//		The final chunk is processed byte-by-byte, so that we do not access
//		any bytes beyond the bit arrays.
//
//----------

void bitwise_fill
   (void*		dstBits,
	const int	bitVal,
	const u64	numBits)
	{
	u64*		dst = (u64*) dstBits;
	u64			n;
	u64			chunk64 = 0;

	if ((bitVal & 1) == 1)
		chunk64 = (u64) (-1);

	for (n=numBits ; n>=64 ; n-=64,dst++)
		*dst = chunk64;

	if (n == 0) return;

	u8*	dstb   = (u8*) dst;
	u8	chunk8 = (u8) chunk64;
	for ( ; n>=8 ; n-=8,dstb++)
		*dstb = chunk8;

	if (n == 0) return;

	u8 mask = least_significant(u8,n);
	*dstb = (chunk8 & mask) | (*dstb & ~mask);
	}

//----------
//
// bitwise_squeeze--
//	Make a copy of a bit array, only copying specified bit positions and
//	packing the result, so that unspecified bit positions are squeezed out.
//
// Example:
//	bits:        11000110001110001101000111011100111010010010001010
//	specBits:    01110011011100100110111111000100110010001001100100
//	copied bits: -100--10-011--0--10-000111---1--11--1---0--00--0--
//	result:      1001001101000011111110000
//
//----------
//
// Arguments:
//	const void*	bits:		Bit array to read.
//	const void*	specBits:	Bit array indicating the bit positions to copy.
//							.. 1s indicate positions to be copied.
//							.. 0s indicate positions to be ignored.
//	u64			numBits:	The length of the input bit arrays (bits and
//							.. specBits), counted in *bits*. See note (1) below.
//	void*		dstBits:	Bit array to fill.  Note that this must be as long
//							.. as bits and specBits.
//	u64			numDstBits:	The length of the output bit array (dstBits),
//							.. counted in *bits*. By default we assume this is
//							.. the same as the input bit arrays.
//
// Returns:
//	The number of bits copied to dstBits; note that the rest of dstBits has
//	been cleared.  It is up to the caller to modify whatever data object
//	contains dstBits so that it knows the new length (if that's important).
//	Note that the return value is never greater than numDstBits.
//
//----------
//
// Notes:
//	(1)	The number of bytes in the bit arrays is ceil(numBits/8) or
//		ceil(numDstBits/8). When numBits is not a multiple of 8, the remaining
//		bits are read from the least significant bits of the final byte.
//	(2)	We process the bytes in 64-bit chunks until we get to the final chunk.
//		The final chunk is processed byte-by-byte, so that we do not access
//		any bytes beyond the bit arrays.
//
//----------

std::uint64_t bitwise_squeeze
   (const void*	bits,
	const void*	specBits,
	const u64	numBits,
	void*		dstBits,
	const u64	_numDstBits)
	{
	u64*		src  = (u64*) bits;
	u64*		scan = (u64*) specBits;
	u64*		dst  = (u64*) dstBits;
	u64			numDstBits;
	u64			n;
	u8*			scanb, *srcb, *chunkb, *dstb;
	u64			bitsWritten = 0;  // placate compiler, re goto hopping over true initialization

	if (_numDstBits == ((u64)-1)) numDstBits = numBits;
	                         else numDstBits = _numDstBits;

	// copy from full 64-bit chunks

	u64 dstChunk    = 0;
	u64 bitsInChunk = 0;
	u64 bitsInDst   = 0;
	for (n=numBits ; n>=64 ; n-=64)
		{
		u64 specChunk = *(scan++);
		u64 srcChunk  = *(src++);

		for (u8 b=0 ; b<64 ; b++)
			{
			if ((specChunk&1) == 1)
				{
				dstChunk |= (srcChunk&1) << bitsInChunk;
				if (++bitsInChunk == 64)
					{
					if (bitsInDst+64 > numDstBits) goto overrun;
					*(dst++) = dstChunk;
					dstChunk    = 0;
					bitsInChunk = 0;
					bitsInDst   += 64;
					}
				}
			specChunk >>= 1;
			srcChunk  >>= 1;
			}
		}

	// copy from full bytes, if any remain

	scanb = (u8*) scan;
	srcb  = (u8*) src;

	if (n > 0)
		{
		for ( ; n>=8 ; n-=8)
			{
			u8 specByte = *(scanb++);
			u8 srcByte  = *(srcb++);

			for (u8 b=0 ; b<8 ; b++)
				{
				if ((specByte&1) == 1)
					{
					dstChunk |= ((u64) (srcByte&1)) << bitsInChunk;
					if (++bitsInChunk == 64)
						{
						if (bitsInDst+64 > numDstBits) goto overrun;
						*(dst++) = dstChunk;
						dstChunk    = 0;
						bitsInChunk = 0;
						bitsInDst   += 64;
						}
					}
				specByte >>= 1;
				srcByte  >>= 1;
				}
			}
		}

	// copy from partial byte, if one remains

	if (n > 0)
		{
		u8 specByte = *scanb;
		u8 srcByte  = *srcb;

		for (u8 b=0 ; b<n ; b++)
			{
			if ((specByte&1) == 1)
				{
				dstChunk |= ((u64) (srcByte&1)) << bitsInChunk;
				if (++bitsInChunk == 64)
					{
					if (bitsInDst+64 > numDstBits) goto overrun;
					*(dst++) = dstChunk;
					dstChunk    = 0;
					bitsInChunk = 0;
					bitsInDst   += 64;
					}
				}
			specByte >>= 1;
			srcByte  >>= 1;
			}
		}

	// write partial chunk, if we have one; note that if we're writing into a
	// final partial chunk of the destination vector, we also clear the rest
	// of that partial chunk and return

	if (bitsInDst + bitsInChunk > numDstBits) goto overrun;

	bitsWritten = bitsInDst;

	if (bitsInChunk > 0)
		{
		if (numDstBits-bitsInDst >= 64)
			{
			*(dst++) = dstChunk;
			bitsInDst   += bitsInChunk;
			bitsWritten += 64;
			}
		else
			{
			u64 bytesLeft = (numDstBits+7-bitsInDst) / 8;
			chunkb = (u8*) &dstChunk;
			dstb   = (u8*) dst;

			for (n=bitsInChunk ; n>=8 ; n-=8)
				{ *(dstb++) = *(chunkb++); bytesLeft--; }
			if (n > 0)
				{ *(dstb++) = *(chunkb++); bytesLeft--; }
			while (bytesLeft > 0)  // erase remaining bytes in this partial chunk
				{ *(dstb++) = 0;           bytesLeft--; }
			bitsInDst += bitsInChunk;

			return bitsInDst;
			}
		}

	// erase full 64-bit chunks, if any remain

	for ( ; bitsWritten+64<numDstBits ; bitsWritten+=64)
		{ *(dst++) = 0;  bitsWritten += 64; }

	// erase bytes, if any remain

	dstb = (u8*) dst;
	for ( ; bitsWritten<numDstBits ; bitsWritten+=8)
		{ *(dstb++) = 0;  bitsWritten += 8; }

	return bitsInDst;

	// we have a full or partial 64-bit chunk but we don't have enough bits
	// left in the destination vector;  write some of the chunk, byte-by-byte

//……… need to test this

overrun:
	if (bitsInDst == numDstBits)
		return bitsInDst;

	u64 bitsToWrite = numDstBits - bitsInDst;   // guaranteed to be less than 64
	if (bitsInChunk > bitsToWrite) bitsInChunk = bitsToWrite;
	dstChunk &= (((u64)1) << bitsToWrite) - 1;  // erase excess bits from ms end of chunk

	chunkb = (u8*) &dstChunk;
	dstb   = (u8*) dst;
	for ( ; bitsToWrite>=8 ; bitsToWrite-=8)
		*(dstb++) = *(chunkb++);
	if (bitsToWrite > 0)
		*(dstb++) = *(chunkb++);

	return bitsInDst + bitsInChunk;
	}

//----------
//
// bitwise_count --
//	Compute the count of 1s in a bit array.
//
//----------
//
// Arguments:
//	const void*	bits:		Bit array to read.
//	u64			numBits:	The length of the bit array, counted in *bits*.
//							.. See note (1) below.
//
// Returns:
//	The number of positions that are 1 in the bit array.
//
//----------
//
// Notes:
//	(1)	The number of bytes in the bit arrays is ceil(numBits/8). When numBits
//		is not a multiple of 8, the remaining bits are read from the least
//		significant bits of the final byte.
//
//----------

u64 bitwise_count
   (const void*	bits,
	const u64	numBits)
	{
	u8*			scan = (u8*) bits;
	u64			n;
	u64			numOnes = 0;

	for (n=numBits ; n>=8 ; n-=8)
		numOnes += popCount8[*(scan++)];

	if (n == 0) return numOnes;

	u8 mask = least_significant(u8,n);
	numOnes += popCount8[(*scan) & mask];

	return numOnes;
	}

//----------
//
// bitwise_and_count, bitwise_mask_count, bitwise_or_count,
// bitwise_or_not_count, and bitwise_xor_count --
//	Compute the count of 1s in the result of a bitwise operation between two
//	bit arrays.
//
//----------
//
// Arguments:
//	const void*	bits1, bits2:	Bit arrays to compare.
//	u64			numBits:		The length of the bit arrays, counted in *bits*.
//								.. See note (1) below.
//
// Returns:
//	The number of positions that are 1 in the result of the bitwise operation.
//
//----------
//
// Notes:
//	(1)	The number of bytes in the bit arrays is ceil(numBits/8). When numBits
//		is not a multiple of 8, the remaining bits are read from the least
//		significant bits of the final byte.
//	(2)	bitwise_xor_count() is equivalent to hamming distance.
//
//----------

u64 bitwise_and_count
   (const void*	bits1,
	const void*	bits2,
	const u64	numBits)
	{
	u8*			scan1 = (u8*) bits1;
	u8*			scan2 = (u8*) bits2;
	u64			n;
	u64			numOnes = 0;

	for (n=numBits ; n>=8 ; n-=8)
		numOnes += popCount8[*(scan1++) & *(scan2++)];

	if (n == 0) return numOnes;

	u8 mask = least_significant(u8,n);
	numOnes += popCount8[(*scan1 & *scan2) & mask];

	return numOnes;
	}

u64 bitwise_mask_count
   (const void*	bits1,
	const void*	bits2,
	const u64	numBits)
	{
	u8*			scan1 = (u8*) bits1;
	u8*			scan2 = (u8*) bits2;
	u64			n;
	u64			numOnes = 0;

	for (n=numBits ; n>=8 ; n-=8)
		numOnes += popCount8[*(scan1++) & ~*(scan2++)];

	if (n == 0) return numOnes;

	u8 mask = least_significant(u8,n);
	numOnes += popCount8[(*scan1 & ~*scan2) & mask];

	return numOnes;
	}

u64 bitwise_or_count
   (const void*	bits1,
	const void*	bits2,
	const u64	numBits)
	{
	u8*			scan1 = (u8*) bits1;
	u8*			scan2 = (u8*) bits2;
	u64			n;
	u64			numOnes = 0;

	for (n=numBits ; n>=8 ; n-=8)
		numOnes += popCount8[*(scan1++) | *(scan2++)];

	if (n == 0) return numOnes;

	u8 mask = least_significant(u8,n);
	numOnes += popCount8[(*scan1 | *scan2) & mask];

	return numOnes;
	}

u64 bitwise_or_not_count
   (const void*	bits1,
	const void*	bits2,
	const u64	numBits)
	{
	u8*			scan1 = (u8*) bits1;
	u8*			scan2 = (u8*) bits2;
	u64			n;
	u64			numOnes = 0;

	for (n=numBits ; n>=8 ; n-=8)
		numOnes += popCount8[*(scan1++) | ~*(scan2++)];

	if (n == 0) return numOnes;

	u8 mask = least_significant(u8,n);
	numOnes += popCount8[(*scan1 | ~*scan2) & mask];

	return numOnes;
	}

u64 bitwise_xor_count
   (const void*	bits1,
	const void*	bits2,
	const u64	numBits)
	{
	u8*			scan1 = (u8*) bits1;
	u8*			scan2 = (u8*) bits2;
	u64			n;
	u64			numOnes = 0;

	for (n=numBits ; n>=8 ; n-=8)
		numOnes += popCount8[*(scan1++) ^ *(scan2++)];

	if (n == 0) return numOnes;

	u8 mask = least_significant(u8,n);
	numOnes += popCount8[(*scan1 ^ *scan2) & mask];

	return numOnes;
	}
