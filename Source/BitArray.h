
#pragma once


#include <sdla/Exception.h>
#include <cstring>


// Simple class that manages arrays of bits, packing them into an allocated
// unsigned char array.
class cBitArray
{
public:
	typedef unsigned char	Byte;

	class Modifier
	{
	public:
		Modifier(cBitArray& array, const int index) :

			m_BitArray(array),
			m_Index(index)

		{
		}


		// Set bit to true or false
		Modifier& operator = (const int value)
		{
			// Calculate bit
			int bit = 1 << (m_Index & 7);

			// First clear the entry
			m_BitArray.GetByte(m_Index >> 3) &= ~bit;

			// OR the value into the array
			m_BitArray.GetByte(m_Index >> 3) |= bit * value;

			return (*this);
		}


		// Convert bit state to integer
		operator int (void) const
		{
			return (!!(m_BitArray.GetByte(m_Index >> 3) & (1 << (m_Index & 7))));
		}

	private:
		// Array to modify
		cBitArray&	m_BitArray;

		// Bit index
		int	m_Index;
	};

	// Can only construct from the bit count
	cBitArray(const int nb_bits) :

		m_Data(0),
		m_NbBytes((nb_bits + 7) >> 3),
		m_NbBits(nb_bits)

	{
		// Allocate enough space for all the bits
		m_Data = new Byte[m_NbBytes];
		memset(m_Data, 0, m_NbBytes);
	}


	~cBitArray(void)
	{
		if (m_Data)
			delete [] m_Data;
	}


	// Get byte for read
	Byte GetByte(const int index) const
	{
		ASSERT(index >= 0 && index < m_NbBytes);
		return (m_Data[index]);
	}


	// Get byte for modification
	Byte& GetByte(const int index)
	{
		ASSERT(index >= 0 && index < m_NbBytes);
		return (m_Data[index]);
	}


	// Check to see if a bit is set
	int operator [] (const int index) const
	{
		ASSERT(index >= 0 && index < m_NbBits);
		return (!!(m_Data[index >> 3] & (1 << (index & 7))));
	}


	// Get a bit for modifying
	Modifier operator [] (const int index)
	{
		ASSERT(index >= 0 && index < m_NbBits);
		return (Modifier(*this, index));
	}


	int GetNbBits(void) const
	{
		return (m_NbBits);
	}


	int GetNbBytes(void) const
	{
		return (m_NbBytes);
	}

private:
	// Allocated array of bits
	Byte*	m_Data;

	// Size of the array
	int		m_NbBytes;
	int		m_NbBits;
};