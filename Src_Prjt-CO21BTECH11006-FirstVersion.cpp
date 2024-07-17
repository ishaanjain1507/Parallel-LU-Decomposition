#include <bits/stdc++.h>

using namespace std::chrono;
using namespace std;

vector<vector<double>> matrix;
int n_procs, mat_size, r; // r is no. of columns that each processor stores locally

vector<double> ReceiveWaitTimes;
atomic<int> ReceiveCount;

// FUNCTION DECLARATIONS
void readInput();
void runner(int id);
int alloc(int k);

class ProcessorNode
{
public:
    int proc_id;
    ProcessorNode *next;

    vector<double> buffer;
    vector<vector<double>> localCols; // local columns that each thread have
    volatile atomic<bool> received;   // variable used to synchronize communication b/w threads
    int l;                            // local index for columns
    int r_local;

    ProcessorNode(int id, int n)
    {
        proc_id = id;
        next = NULL;
        buffer.resize(n);
        localCols.resize(n);
        l = 0;
        received.store(false);
        r_local = 0;
    }

    void Send()
    {
        cout << "In the send " << proc_id << endl;
        for (int i = 0; i < buffer.size(); i++)
        {
            next->buffer[i] = buffer[i];
        }
        this->next->received.store(true);
        cout << "Going out of send " << proc_id << endl;
        return;
    }

    void Receive()
    {
        cout << "In the receive " << proc_id << endl;

        auto startTime = high_resolution_clock::now();

        while (!this->received.load())
        {
        } // spin until it receives the buffer from previous thread

        this->received.store(false);

        auto endTime = high_resolution_clock::now();
        auto TotalTime = (int)duration_cast<microseconds>(endTime - startTime).count();

        ReceiveWaitTimes[proc_id] += TotalTime;
        ReceiveCount.fetch_add(1);

        for (int i = 0; i < buffer.size(); i++)
            cout << buffer[i] << endl;
        cout << "Going out of the receive " << proc_id << endl;
        return;
    }

    void Prep(int k)
    {
        cout << "In prep: Proc=" << proc_id << ", k=" << k << endl;
        cout << "    Dividing by " << localCols[k][l] << endl;

        for (int i = k + 1; i < mat_size; i++)
        {
            localCols[i][l] = -localCols[i][l] / localCols[k][l];
        }

        for (int i = 0; i < mat_size; i++)
            this->buffer[i] = localCols[i][l];
        cout << "Going out of prep: Proc=" << proc_id << ", k=" << k << endl;
        return;
    }

    void Update(int k, int j)
    {
        cout << "In update: Proc=" << proc_id << ", j=" << j << ", k=" << k << endl;

        for (int i = k + 1; i < mat_size; i++)
        {
            localCols[i][j] = localCols[i][j] + buffer[i] * localCols[k][j];
        }

        cout << "Going out of update: Proc=" << proc_id << ", j=" << j << ", k=" << k << endl;
        return;
    }

    void BROADCAST(int col_ind)
    {
        int k = alloc(col_ind);
        int pred_id = (k != 0) ? k - 1 : n_procs - 1;
        if (proc_id == k)
        {
            Send();
        }
        else
        {
            if (proc_id == pred_id)
            {
                cout << "In broadcast second if " << proc_id << endl;
                Receive();
            }
            else
            {
                cout << "In broadcast last else " << proc_id << endl;
                Receive();
                Send();
            }
        }
    }

    void FACTORIZATION()
    {
        for (int k = 0; k < mat_size - 1; k++)
        {
            if (alloc(k) == proc_id)
            {
                Prep(k);
                l++;
                Send();
            }
            // BROADCAST(alloc(k), buffer);

            else
            {
                Receive();
            }
            for (int j = l; j < r_local; j++)
            {
                Update(k, j);
            }
            if (alloc(k) != proc_id && proc_id != (k - 1) % n_procs)
            {
                Send();
            }
        }
    }

    void printNode()
    {
        cout << "  Processor ID: " << proc_id << endl;
        for (auto it : localCols)
        {
            cout << "    ";
            for (auto pt : it)
                cout << pt << " ";
            cout << endl;
        }
        return;
    }
};

class CircularLinkedList
{
public:
    int numprocs;
    ProcessorNode *head;
    vector<ProcessorNode *> procs;
    vector<vector<double>> finalMatrix;

    CircularLinkedList(int numprocs, int matsize)
    {
        this->numprocs = numprocs;
        ProcessorNode *prev = new ProcessorNode(0, matsize);
        procs.push_back(prev);
        this->head = prev;
        for (int i = 1; i < numprocs; i++)
        {
            ProcessorNode *node = new ProcessorNode(i, matsize);
            procs.push_back(node);
            prev->next = node;
            prev = node;
        }
        prev->next = this->head;
    }

    void DivideAmongThreads()
    {
        ProcessorNode *temp;

        // storing columns locally in each proccessor node
        for (int col = 0; col < mat_size; col++)
        {
            int k = alloc(col);
            temp = procs[k];
            for (int row = 0; row < mat_size; row++)
            {
                temp->localCols[row].push_back(matrix[row][col]);
            }
            temp->r_local++;
        }
        return;
    }

    void verifyAllocation()
    {
        cout << "verifyAllocation(): " << endl;
        ProcessorNode *curr = head;
        do
        {
            curr->printNode();
            curr = curr->next;
        } while (curr != head);
        return;
    }

    void generateOutput()
    {
        finalMatrix.resize(mat_size, vector<double>(mat_size));
        for (int j = 0; j < mat_size; j++)
        {
            int k = alloc(j);
            ProcessorNode *temp = procs[k];
            temp->l = 0;
            for (int i = 0; i < mat_size; i++)
            {
                finalMatrix[i][j] = temp->localCols[i][temp->l];
            }
            temp->l++;
        }
        fstream output;
        output.open("LU-FirstVersion-Output.txt", ios::out);

        if (!output)
        {
            cout << "File couldn't be opened\n";
        }
        output << "Final matrix: " << endl;
        for (int i = 0; i < mat_size; i++)
        {
            output << "  ";
            for (int j = 0; j < mat_size; j++)
            {
                output << finalMatrix[i][j] << " ";
            }
            output << endl;
        }
        output << endl;
        output << "Lower triangular matrix: " << endl;
        for (int i = 0; i < mat_size; i++)
        {
            output << "  ";
            for (int j = 0; j < mat_size; j++)
            {
                if (i == j)
                    output << "1 ";
                else if (j > i)
                    output << "0 ";
                else
                    output << -finalMatrix[i][j] << " ";
            }
            output << endl;
        }
        output << endl;
        output << "Upper triangular matrix: " << endl;
        for (int i = 0; i < mat_size; i++)
        {
            output << "  ";
            for (int j = 0; j < mat_size; j++)
            {
                if (j >= i)
                    output << finalMatrix[i][j] << " ";
                else
                    output << "0 ";
            }
            output << endl;
        }
        output << endl;
        output.close();
    }
};

CircularLinkedList *CircularNodes;
pthread_barrier_t b;

// MAIN FUNCTION
int main()
{

    // Step-1: Read the input matrix
    readInput();
    cout << "Number of processors: " << n_procs << endl;
    cout << "Matrix size: " << mat_size << endl;
    cout << "Input matrix: " << endl;
    for (auto it : matrix)
    {
        cout << "  ";
        for (auto pt : it)
            cout << pt << " ";
        cout << endl;
    }

    r = (int)(mat_size / n_procs);

    // Step-2: Create n_procs
    CircularNodes = new CircularLinkedList(n_procs, mat_size);

    // Initialize barrier
    pthread_barrier_init(&b, NULL, n_procs);

    // For Receive wait time
    ReceiveWaitTimes.resize(n_procs);
    ReceiveCount.store(0);

    // Step-3: Allocate cols to each of the processors
    CircularNodes->DivideAmongThreads();
    // Testing
    CircularNodes->verifyAllocation();

    auto startTime = high_resolution_clock::now();
    // Step-4: Create 'n_procs' number of threads
    vector<thread> proc_threads;
    for (int i = 0; i < n_procs; i++)
    {
        proc_threads.push_back(thread(runner, i));
    }

    // Step-5: Wait till all threads terminate
    for (auto &thr : proc_threads)
    {
        thr.join();
    }

    // Destory barrier
    pthread_barrier_destroy(&b);

    cout << "\nAll threads terminated.\nFinal matrix:" << endl;

    auto endTime = high_resolution_clock::now();

    // Step-6: Write to output file
    CircularNodes->generateOutput();
    CircularNodes->verifyAllocation();

    auto TotalTime = (int)duration_cast<microseconds>(endTime - startTime).count();
    cout << "Total time taken: " << TotalTime << "ms" << endl;

    cout << "Throughput: " << (double)n_procs / (double)TotalTime << endl;

    double AvgWaitTime = (double)accumulate(ReceiveWaitTimes.begin(), ReceiveWaitTimes.end(), 0) / (double)ReceiveCount.load();
    cout << "Average Wait Time in Receive: " << AvgWaitTime << "ms" << endl;

    return 0;
}

void runner(int id)
{
    ProcessorNode *myProc = CircularNodes->procs[id];

    pthread_barrier_wait(&b); // wait until all threads get their processor nodes

    // factorization
    myProc->FACTORIZATION();
    return;
}

// readInput reads the input file "input.txt"
// to get values n_procs, mat_size and matrix (global vars)
void readInput()
{
    ifstream input("inp-params.txt");
    if (!input)
    {
        cout << "readInput(): Input file could not be opened." << endl;
        return;
    }
    input >> n_procs >> mat_size;
    int n = mat_size;
    matrix.resize(n);
    for (int i = 0; i < n; i++)
    {
        matrix[i].resize(n);
        for (int j = 0; j < n; j++)
            input >> matrix[i][j];
    }
    cout << "readInput(): Input file read and data saved." << endl;
    input.close();
    return;
}

int alloc(int k)
{
    return k % n_procs;
}