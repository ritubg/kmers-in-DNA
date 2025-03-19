# kmers-in-DNA
This is a C++ project mainly focused on finding k-mers in DNA sequence using Hash table data structure

---

## 🧬 Overview
This is a C++ implementation of a **hash table** designed for efficient processing of **DNA k-mers**.  
A **k-mer** is a substring of length `k` from a DNA sequence, which consists of four nucleotides:  

- **Adenine (A)**  
- **Thymine (T)**  
- **Guanine (G)**  
- **Cytosine (C)**  

The **DNA sequence** is stored in a file named `dna.txt`.  
This program efficiently indexes the genome using a hash table and provides various operations.

---

## 🚀 Features
✔ **Create a hash table** for storing k-mers.  
✔ **Insert and delete** k-mers using **linear probing** for collision resolution.  
✔ **Replace** a k-mer in the hash table.  
✔ **Find the location** of a specific k-mer.  
✔ **Calculate the frequency** of a specific k-mer.  
✔ **Display** the entire hash table.  
✔ **Automatic resizing** when the occupancy ratio exceeds 50% (table size doubles).  
✔ **Exception handling** for invalid k-mer sizes and out-of-bounds errors.  

---

## 🔢 Hash Function
The hash function computes the **hash value** by summing the **ASCII values** of the characters in the k-mer  
and taking the **modulus with the table size**:



📌 **Initial table size** = `100`.  

---

## ⏳ Time Complexity  

| Operation           | Time Complexity |
|---------------------|----------------|
| **Hash function**      | O(k) |
| **Search**            | O(k) (linear probing) |
| **Frequency Calculation** | O(k) |
| **Insertion**         | O(k) |
| **Rehashing**        | O(nk) |
| **Deletion**         | O(n) |
| **Display**          | O(n) |
| **Modification (replace k-mers)** | O(nk) |  

---

## 🛠 Installation & Usage  

### 1️⃣ Clone the Repository  
```sh
git clone https://github.com/your-username/kmers-in-DNA.git
cd kmers-in-DNA

### 2️⃣ Compile the Program
g++ -o kmerhash "main(2).cpp"


### 3️⃣ Run the Program
./kmerhash

Ensure dna.txt is in the same directory as the executable.

