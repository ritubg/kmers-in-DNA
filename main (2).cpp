#include <iostream>
#include <vector>
#include<string>
#include<fstream>
#include <algorithm>

using namespace std;

const int MAXLOC = 150;
class KmerHash
{
    private:
        int tableSize;
        string *Hash_table;
        int**locations;
        int*locationCount;
        int numOfUniqueKeys;
    public:
        KmerHash(int m)
        {
            tableSize = m;
            numOfUniqueKeys = 0;
            Hash_table = new string[tableSize];
            locations = new int*[tableSize];
            for(int i=0;i<tableSize;i++)
                locations[i] = new int[MAXLOC];
            locationCount = new int[tableSize]();
        }
        int Hash_function(const string& key)
        {
            int hash=0;
            for(int i=0;i<key.length();i++)
            {
                hash = hash + key[i];
            }
            int index=hash%tableSize;
            return index;
        }
        void insertKmer(const string&kmer,int pos)
        {
            int index=Hash_function(kmer);
            while (!Hash_table[index].empty()&& Hash_table[index] != kmer)
            {
                index = (index + 1) % tableSize;
            }

            if(Hash_table[index] == kmer)
            {
                if(locationCount[index] < MAXLOC)
                {
                    int x=locationCount[index];
                    locations[index][x] = pos;
                    locationCount[index]++;
                }
                else
                {
                    cout<<"MAXIMUM LOCATIONS REACHED FOR Kmer "<<kmer<<endl;
                }
            }
            else
            {
                Hash_table[index]=kmer;
                locations[index][0] = pos;
                locationCount[index]=1;
                numOfUniqueKeys++;
            }

            float occupancy = (float)numOfUniqueKeys/tableSize;
            if(occupancy>0.5)
            {
                rehashData();
            }
        }

        void rehashData()
        {
            int oldTableSize = tableSize;
            string *oldHashTable =  Hash_table;
            int**oldLocations = locations;
            int*oldLocationCount = locationCount;

            tableSize = tableSize*2;
            numOfUniqueKeys = 0;
            Hash_table = new string[tableSize];
            locations = new int*[tableSize];
            for(int i=0;i<tableSize;i++)
                locations[i] = new int[MAXLOC];
            locationCount = new int[tableSize]();


            for (int i = 0; i < oldTableSize; i++)
            {
                if (!oldHashTable[i].empty())
                {
                    for (int j = 0; j < oldLocationCount[i]; j++)
                    {
                    insertKmer(oldHashTable[i], oldLocations[i][j]);
                    }
                }
            }
            for(int i=0;i<oldTableSize;i++)
                delete[] oldLocations[i];
            delete[] oldLocations;
            delete[] oldHashTable;
            delete[] oldLocationCount;

        }


        ~KmerHash()
        {
            for(int i=0;i<tableSize;i++)
                delete[] locations[i];
            delete[] locations;
            delete[] Hash_table;
            delete[] locationCount;
        }

        vector<int>locationOfSpecificKmer(const string&kmer)
        {
            vector<int>locs;
            int index=Hash_function(kmer);
            while (!Hash_table[index].empty() && Hash_table[index] != kmer)
            {
                index = (index + 1) % tableSize;
            }

            if(Hash_table[index] == kmer)
            {
                for(int i=0;i<locationCount[index];i++)
                {
                    locs.push_back(locations[index][i]);
                }
            }
            return locs;
        }


        void printKmers()
        {
            cout<<"\t\t\tK-mers"<<endl;
            for(int i=0;i<tableSize;i++)
            {
                if(!Hash_table[i].empty())
                {
                    cout<<Hash_table[i]<<endl;
                }
            }
        }

        void displayHashTable()
        {
            cout<<"\n\nK-mer DNA Hash table is\n\n";
            for(int i=0;i<tableSize;i++)
            {
                if(!Hash_table[i].empty())
                {
                    cout<<" | "<<i<<"\t| "<<"K-mer"<<" "<<Hash_table[i]<<" | "<<" LOCATIONS : ";
                    for(int j=0;j<locationCount[i];j++)
                    {
                        cout<<locations[i][j]<<" ";
                    }
                    cout<<endl;
                }
            }
        }
        void deletekmer(const string &kmer, string &genome, int k)
        {
            int index = Hash_function(kmer);
            bool found = false;

            while (!Hash_table[index].empty() && Hash_table[index] != kmer)
            {
                index = (index + 1) % tableSize;
            }

            if (Hash_table[index] == kmer)
            {
                found = true;
                int count = locationCount[index];
                vector<int> positions(locations[index], locations[index] + count);
                sort(positions.rbegin(), positions.rend());

                for (int pos : positions)
                {
                    if (pos + kmer.length() <= genome.length())
                    {
                        genome.erase(pos, kmer.length());
                    }
                }

                Hash_table[index].clear();
                locationCount[index] = 0;
                for (int i = 0; i < MAXLOC; i++)
                {
                    locations[index][i] = 0;
                }

                cout << "\nSuccessfully deleted all occurrences of the k-mer: " << kmer << endl;
                cout << "\n----------------------------------------- ********* -----------------------------------------" << endl;
            }

            if (found)
            {
                ofstream fileo("dna.txt");
                fileo << genome;
                fileo.close();
            }
            else
            {
                cout << "\nThe given k-mer: " << kmer << " not found " << endl;
                cout << "\n----------------------------------------- ********* -----------------------------------------" << endl;
            }

            string *oldHashTable = Hash_table;
            int **oldLocations = locations;
            int *oldLocationCount = locationCount;

            numOfUniqueKeys = 0;
            Hash_table = new string[tableSize];
            locations = new int *[tableSize];
            for (int i = 0; i < tableSize; i++)
            locations[i] = new int[MAXLOC];
            locationCount = new int[tableSize]();

            ifstream file("dna.txt");
            if (!file.is_open())
            {
                cout << "\nError opening file!" << endl;
                exit(0);
            }

            string genomeNew, line;
            while (getline(file, line))
            {
                genomeNew = genomeNew + line;
            }
            file.close();

            for (int i = 0; i <= genomeNew.length() - k; ++i)
            {
                string kmer = genomeNew.substr(i, k);
                insertKmer(kmer, i);
            }

            for (int i = 0; i < tableSize; i++)
            delete[] oldLocations[i];
            delete[] oldLocations;
            delete[] oldHashTable;
            delete[] oldLocationCount;
        }


        int frequencyOfkmer(const string& kmer)
        {
            int index=Hash_function(kmer);
            while(!Hash_table[index].empty() && Hash_table[index]!=kmer)
            {
                index=(index+1)%tableSize;
            }
            if(Hash_table[index]==kmer)
            {
                return locationCount[index];
            }
            return 0;
        }
        bool searchKmer(const string& kmer)
        {
            int index= Hash_function(kmer);
            while(!Hash_table[index].empty() && Hash_table[index]!=kmer)
            {
                index=(index+1)%tableSize;
            }
            if(Hash_table[index]==kmer)
            {
                return true;
            }
            return false;
        }
        void replaceKmer(const string &oldKmer,const string &newKmer,string &genome,int k)
        {
            int index=Hash_function(oldKmer);
            bool found=false;
            while(!Hash_table[index].empty() && Hash_table[index]!=oldKmer)
            {
                index=(index+1)%tableSize;
            }
            if(Hash_table[index]==oldKmer)
            {
                found=true;
                int count=locationCount[index];
                for(int i=0;i<count;i++)
                {
                    int pos=locations[index][i];
                    genome.replace(pos,oldKmer.length(),newKmer);
                }
                Hash_table[index]=newKmer;
                cout<<"\nSuccessfully replaced the K-mer "<<oldKmer<<" with "<<newKmer<<endl;
            }
            if(found)
            {
                ofstream fileo("dna.txt");
                fileo<<genome;
                fileo.close();
            }
            else
            {
                cout<<"\nK-mer "<<oldKmer<<" not found in the hash table."<<endl;
            }

            string *oldHashTable =  Hash_table;
            int**oldLocations = locations;
            int*oldLocationCount = locationCount;

            numOfUniqueKeys = 0;
            Hash_table = new string[tableSize];
            locations = new int*[tableSize];
            for(int i=0;i<tableSize;i++)
                locations[i] = new int[MAXLOC];
            locationCount = new int[tableSize]();


            ifstream file("dna.txt");
            if (!file.is_open())
            {
                cout << "\nError opening file!" << endl;
                exit(0);
            }

            string genomeNew,line;
            while (getline(file, line))
            {
                genomeNew = genomeNew + line;
            }
            file.close();

            for (int i = 0; i <= genomeNew.length() - k; ++i)
            {
                string kmer = genomeNew.substr(i, k);
                insertKmer(kmer, i);
            }


            for(int i=0;i<tableSize;i++)
                delete[] oldLocations[i];
            delete[] oldLocations;
            delete[] oldHashTable;
            delete[] oldLocationCount;
        }

};


bool isValidKmer(const string &dna)
{
    for(char nucleotide : dna)
    {
        if(nucleotide!='A' && nucleotide!='T' && nucleotide!='C' && nucleotide!='G')
        {
            return false;
        }
    }
    return true;
}

int main()
{
    int choice1,choice2,choice3;
    cout<<"\t\t\t\t\t----------------------------------------------"<<endl;
    cout<<"\t\t\t\t\t                      MENU                    "<<endl;
    cout<<"\t\t\t\t\t----------------------------------------------"<<endl;
    cout<<"\n1. Basic operations in File.\n\n2. K-mers operations in Hash-Table.\n";
    cout<<"\nEnter your choice\t: ";
    cin>>choice1;
    cout<<"\n-----------------------------------------             *********             -----------------------------------------"<<endl;
    if(choice1==1)
    {
        while(1)
        {
            ifstream file("dna.txt");
            if (!file.is_open())
            {
                cout << "Error in opening a file!" << endl;
                return 1;
            }

            string genome,line;
            while (getline(file, line))
            {
                genome = genome + line;
            }
            file.close();
            cout<<"\n\n1. Insert DNA genome to a specific location in file\n\n2. Delete DNA genome for a specific length in file\n\n3. Display the file contents";
            cout<<"\n\nEnter your choice\t: ";
            cin>>choice2;
            if(choice2 == 1)
            {
                int index;
                string dna;
                cout<<"\nEnter DNA genome : ";
                cin>>dna;
                if(!isValidKmer(dna))
                {
                    cout<<"\nInvalid input...The K-mer must contain only A,T,G or C"<<endl;
                    cout<<"\n-----------------------------------------             *********             -----------------------------------------"<<endl;
                }
                else
                {
                    cout<<"\nEnter location to insert : ";
                    cin>>index;
                    genome.insert(index,dna);
                    cout<<"\n"<<dna<<" has been inserted successfully\n";
                    cout<<"\n-----------------------------------------             *********             -----------------------------------------"<<endl;
                    ofstream fileo("dna.txt");
                    fileo<<genome;
                    fileo.close();
                }
            }
            else if(choice2 == 2)
            {
                int options;
                cout<<"\n\n1. Delete at a specific location\n";
                cout<<"\n2. Delete a substring DNA in the genome\n";
                cout<<"\n\nChoose the deleting option \t:";
                cin>>options;
                if(options==1)
                {
                    int index,length;
                    cout<<"\nEnter the start index\t:";
                    cin>>index;
                    cout<<"\nEnter the length of the DNA genome to delete\t:";
                    cin>>length;
                    genome.erase(index,length);
                    ofstream fileo("dna.txt");
                    fileo<<genome;
                    fileo.close();
                    cout<<"\n\nSuccessfully Deleted at the given Index "<<index<<" of length "<<length<<endl;
                    cout<<"\n-----------------------------------------             *********             -----------------------------------------"<<endl;
                }
                else if(options==2)
                {
                    string dna;
                    int index;
                    char option_2;
                    cout << "\nEnter the DNA substring to delete\t: ";
                    cin >> dna;
                    if(!isValidKmer(dna))
                    {
                        cout<<"\nInvalid input...The K-mer must contain only A,T,G or C"<<endl;
                        cout<<"\n-----------------------------------------             *********             -----------------------------------------"<<endl;
                    }
                    else
                    {
                        cout << "\nDo you want to delete a specific occurrence or all occurrences (specific: 'y')\t: ";
                        cin >> option_2;
                        if (option_2 == 'y')
                        {
                            int specific_index;
                            cout << "\n\nEnter the specific index of " << dna << " to delete\t: ";
                            cin >> specific_index;
                            size_t pos = genome.find(dna);
                            int count = 0;
                            while (pos != string::npos)
                            {
                                if (count == specific_index)
                                {
                                    genome.erase(pos, dna.length());
                                    ofstream fileo("dna.txt");
                                    fileo<<genome;
                                    fileo.close();
                                    break;
                                }
                                count++;
                                pos = genome.find(dna, pos + 1);
                            }
                            if (count < specific_index)
                            {
                                cout << "\nThe specific occurrence was not found!" << endl;
                            }
                        }
                        else
                        {
                            size_t pos = genome.find(dna);
                            int count = 0;
                            while (pos != string::npos)
                            {
                                genome.erase(pos, dna.length());
                                count++;
                                pos = genome.find(dna, pos);
                            }
                            cout << "\nDeleted " << count << " occurrences of the substring '" << dna << "' from the genome." << endl;
                            cout<<"\n-----------------------------------------             *********             -----------------------------------------"<<endl;
                            ofstream fileo("dna.txt");
                            fileo<<genome;
                            fileo.close();
                        }
                    }
                }
                else
                {
                    cout<<"\nIn-valid choice"<<endl;
                    cout<<"\n-----------------------------------------             *********             -----------------------------------------"<<endl;
                }
            }
            else if(choice2 == 3)
            {
                cout<<"\nYour File's content is \t:\n"<<endl;
                cout<<genome;
                cout<<"\n-----------------------------------------             *********             -----------------------------------------"<<endl;
                cout<<endl;
            }
            else
            {
                cout<<"\nIn - Valid input"<<endl;
            }
            char ch;
            cout<<"\nWould you like to perform any other operation(y/n)? ";
            cin>>ch;
            if(ch!='y')
            {
                break;
            }
        }
    }
    else if(choice1 == 2)
    {
        ifstream file("dna.txt");
        if (!file.is_open())
        {
            cout << "\nError opening file!" << endl;
            return 1;
        }

        string genome,line;
        while (getline(file, line))
        {
            genome = genome + line;
        }
        file.close();

        int k,genome_length;
        genome_length=genome.length();
        KmerHash table(100);
        cout<<"Enter k value : ";
        cin>>k;
        if(k > genome_length)
        {
            cout<<"\nSorry, "<<k<<" is greater than the length of the genome "<<genome_length<<endl;
            return 1;
        }
        for (int i = 0; i <= genome.length() - k; ++i)
        {
            string kmer = genome.substr(i, k);
            table.insertKmer(kmer, i);
        }
        while(1)
        {
            cout<<"\n1. Print k-mers\n\n2. Find initial locations of a given k-mer using hash table\n\n3. Find frequency of a given k-mer\n";
            cout<<"\n4. Delete a k-mer from all the locations in file\n\n5. Search a kmer in the Hash table \n\n6. Display the k-mer hash table.\n\n7. Modify the existing K-mer";
            cout<<"\n\nEnter your choice : ";
            cin>>choice3;
            if(choice3 == 1)
            {
                cout<<"\nYour K-mers as hash table keys are \t:\n";
                table.printKmers();
                cout<<"\n-----------------------------------------             *********             -----------------------------------------"<<endl;
            }
            else if(choice3 == 2)
            {
                string dna;
                cout<<"\nEnter k-mer : ";
                cin>>dna;
                if(!isValidKmer(dna))
                {
                    cout<<"\nInvalid input...The K-mer must contain only A,T,G or C"<<endl;
                    cout<<"\n-----------------------------------------             *********             -----------------------------------------"<<endl;
                }
                else
                {
                    vector<int>locs = table.locationOfSpecificKmer(dna);
                    if (!locs.empty())
                    {
                        cout<<"\nInitial Locations from the genome file are : ";
                        for (int l : locs)
                        {
                            cout << l << " ";
                        }
                        cout<<"\n-----------------------------------------             *********             -----------------------------------------"<<endl;
                    }
                    else
                    {
                        cout << "\nK-mer "<<dna<<" not found." << endl;
                        cout<<"\n-----------------------------------------             *********             -----------------------------------------"<<endl;
                    }
                }
            }
            else if(choice3 == 3)
            {
                int frequency;
                string dna;
                cout<<"\n\nEnter the DNA of length "<<k<<" to find it's frequency\t:";
                cin>>dna;
                if(!isValidKmer(dna))
                {
                    cout<<"\nInvalid input...The K-mer must contain only A,T,G or C"<<endl;
                    cout<<"\n-----------------------------------------             *********             -----------------------------------------"<<endl;
                }
                else
                {
                    if(dna.length()!=k)
                    {
                        cout<<"\n\nKmer's size is not equal to the size of k"<<endl;
                        cout<<"\n-----------------------------------------             *********             -----------------------------------------"<<endl;
                    }
                    else
                    {
                        frequency=table.frequencyOfkmer(dna);
                        if(frequency==0)
                        {
                            cout<<"\n\nSorry the given kmer : "<<dna<<" is not found in the Hash table"<<endl;
                            cout<<"\n-----------------------------------------             *********             -----------------------------------------"<<endl;

                        }
                        else
                        {
                            cout<<"\nThe Frequency of kmer "<<dna<<" : "<<frequency<<endl;
                            cout<<"\n-----------------------------------------             *********             -----------------------------------------"<<endl;

                        }
                    }
                }
            }
            else if(choice3 == 4)
            {
                string dna;
                cout<<"\nEnter the DNA to delete :";
                cin>>dna;
                if(!isValidKmer(dna))
                {
                    cout<<"\nInvalid input...The K-mer must contain only A,T,G or C"<<endl;
                    cout<<"\n-----------------------------------------             *********             -----------------------------------------"<<endl;
                }
                else
                {
                    if(dna.length()!=k)
                    {
                        cout<<"\nKmer's size is not equal to the size of k"<<endl;
                        cout<<"\n-----------------------------------------             *********             -----------------------------------------"<<endl;

                    }
                    else
                    {
                        table.deletekmer(dna,genome,k);

                    }
                }
            }
            else if(choice3 == 5)
            {
                string dna;
                cout<<"\nEnter the kmer to search in the Genome\t:";
                cin>>dna;
                if(!isValidKmer(dna))
                {
                    cout<<"\nInvalid input...The K-mer must contain only A,T,G or C"<<endl;
                    cout<<"\n-----------------------------------------             *********             -----------------------------------------"<<endl;
                }
                else
                {
                    if(dna.length()!=k)
                    {
                        cout<<"\nKmer's size is not equal to the size of k"<<endl;
                    }
                    else
                    {
                        if(table.searchKmer(dna))
                        {
                            cout<<"\nYES... the given kmer "<<dna<<" is found in the genome"<<endl;
                        }
                        else
                        {
                            cout<<"\nNO... the given kmer "<<dna<<" is not found in the genome"<<endl;
                        }
                    }
                    cout<<"\n-----------------------------------------             *********             -----------------------------------------"<<endl;
                }
            }
            else if(choice3 == 6)
            {
                cout<<"\nYour Hash-Table is \t:\n"<<endl;
                table.displayHashTable();
                cout<<"\n-----------------------------------------             *********             -----------------------------------------"<<endl;

            }
            else if(choice3==7)
            {
                string dna1;
                cout<<"\nEnter the DNA genome to modify\t:";
                cin>>dna1;
                if(!isValidKmer(dna1))
                {
                    cout<<"\nInvalid input...The K-mer must contain only A,T,G or C"<<endl;
                    cout<<"\n-----------------------------------------             *********             -----------------------------------------"<<endl;
                }
                else
                {
                    if(dna1.length()!=k)
                    {
                        cout<<"\nKmer's size is not equal to the size of k"<<endl;
                        cout<<"\n-----------------------------------------             *********             -----------------------------------------"<<endl;
                    }
                    else
                    {
                        string dna2;
                        cout<<"\nEnter the new K-mer\t:";
                        cin>>dna2;
                        if(!isValidKmer(dna2))
                        {
                            cout<<"\nInvalid input...The K-mer must contain only A,T,G or C"<<endl;
                            cout<<"\n-----------------------------------------             *********             -----------------------------------------"<<endl;
                        }
                        else
                        {
                            if(dna2.length()!=k)
                            {
                                cout<<"\nKmer's size is not equal to the size of k"<<endl;
                                cout<<"\n-----------------------------------------             *********             -----------------------------------------"<<endl;
                            }
                            else
                            {
                                table.replaceKmer(dna1,dna2,genome,k);
                            }
                        }
                    }
                }
            }
            else
            {
                cout<<"\nIn - Valid input "<<endl;
            }
            char ch;
            cout<<"\nWould you like to perform any other operation(y/n)? ";
            cin>>ch;
            if(ch!='y')
            {
                 break;
            }
        }
    }
    else
    {
        cout<<"\nIn - Valid input "<<endl;
    }
    return 0;
}

