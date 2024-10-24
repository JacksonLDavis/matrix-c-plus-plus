#include <iostream>

/**
 * Code written by Jackson L. Davis
 *
 * A matrix class for creating basic matrices,
 * only supports matrices of integers
 *
 * Reference Numbers for Exceptions:
 * 501: Matrix() constructor: r must be positive
 * 502: Matrix() constructor: c must be positive
 * 503: getEntry(): 1 <= r <= rows and 1 <= c <= columns must be true
 * 504: setEntry(): 1 <= r <= rows and 1 <= c <= columns must be true
 * 505: getRowMatrix(): 1 <= rowIndex <= rows must be true
 * 506: getColumnMatrix(): 1 <= columnIndex <= columns must be true
 * 507: makeIdentity(): isSquare() must return true
 * 508: determinant(): m.isSquare() must return true
 * 509: sumOfMatrices(): sameDimensions(m1, m2) must return true
 * 510: dotProduct(): canDot(m1, m2) must return true
 * 511: productOfMatrices(): canMultiply(m1, m2) must return true
 */
class Matrix {
    private:
    // int matrix[][]; will not work because C++ needs to know the array dimensions at compile time
    int **matrix; // use a pointer instead
    int rows;
    int columns;

    public:
    /**
     * Constructor for the Matrix class,
     * creates a zero matrix
     *
     * @param r: number of rows for the matrix
     * @param c: number of columns for the matrix
     * @precond: r > 0 && c > 0
     */
    Matrix(int r, int c) {
        // check if r and c are positive
        if (r <= 0) {
            throw 501;
        }
        else if (c <= 0) {
            throw 502;
        }
        // create the Matrix object
        else {
            rows = r;
            columns = c;
            // allocate memory for the array of pointers
            matrix = new int*[rows];
            // allocate memory for each inner array
            for (int i = 0; i < rows; i++) {
                matrix[i] = new int[columns];
                for (int j = 0; j < columns; j++) {
                    matrix[i][j] = 0; // initialize each element to 0
                }
            }
        }
    }

    /**
     * Getter method for rows
     */
    int getRows() {
        return rows;
    }

    /**
     * Getter method for columns
     */
    int getColumns() {
        return columns;
    }

    /**
     * Getter method for matrix
     */
    int** getMatrix() {
        return matrix;
    }

    /**
     * Get the number at a specific entry of the matrix
     * To stay consistent with how matrix entries are notated in mathematics,
     * the top row has an index of 1,
     * and the left-most column has an index of 1
     *
     * @param r: the row of the entry to get
     * @param c: the column of the entry to get
     * @precond: 1 <= r <= rows && 1 <= c <= columns
     * @return: the number at the specified entry
     */
    int getEntry(int r, int c) {
        // subtract 1 from r and c because in programming languages,
        // indexing starts at 0 instead of 1
        r -= 1;
        c -= 1;
        // check for out of bounds indices
        if (r < 0 || c < 0 || r >= rows || c >= columns) {
            throw 503;
        }
        // return the specified entry
        else {
            return matrix[r][c];
        }
    }

    /**
     * Change a specific entry of the matrix
     * To stay consistent with how matrix entries are notated in mathematics,
     * the top row has an index of 1,
     * and the left-most column has an index of 1
     *
     * @param r: the row of the entry to change
     * @param c: the column of the entry to change
     * @param changeTo: the number to change the entry to
     * @precond: 1 <= r <= rows && 1 <= c <= columns
     * @postcond: the matrix entry at row r and column c is changed to changeTo if r and c are in bounds
     */
    void setEntry(int r, int c, int changeTo) {
        // subtract 1 from r and c because in computer science,
        // indexing starts at 0 instead of 1
        r -= 1;
        c -= 1;
        // check for out of bounds indices
        if (r < 0 || c < 0 || r >= rows || c >= columns) {
            throw 504;
        }
        // change the specified entry
        else {
            matrix[r][c] = changeTo;
        }
    }

    /**
     * Get a specific row of a matrix
     * To stay consistent with how matrix entries are notated in mathematics,
     * the top row has an index of 1
     *
     * @param rowIndex: the row index of the matrix to get
     * @precond: 1 <= rowIndex <= rows
     * @return: a 1xn matrix
     */
    Matrix getRowMatrix(int rowIndex) {
        // check for out of bounds index
        if (rowIndex < 1 || rowIndex > rows) {
            throw 505;
        }
        // return the row matrix
        else {
            Matrix rowMatrix(1, columns);
            // i starts at 1 to stay consistent with matrix entry notation
            for (int i = 1; i <= columns; i++) {
                rowMatrix.setEntry(1, i, getEntry(rowIndex, i));
            }
            return rowMatrix;
        }
    }

    /**
     * Get a specific column of a matrix
     * To stay consistent with how matrix entries are notated in mathematics,
     * the left-most column has an index of 1
     *
     * @param columnIndex: the column index of the matrix to get
     * @precond: 1 <= columnIndex <= columns
     * @return: an nx1 matrix
     */
    Matrix getColumnMatrix(int columnIndex) {
        // check for out of bounds index
        if (columnIndex < 1 || columnIndex > columns) {
            throw 506;
        }
        // return the column matrix
        else {
            Matrix columnMatrix(rows, 1);
            // i starts at 1 to stay consistent with matrix entry notation
            for (int i = 1; i <= rows; i++) {
                columnMatrix.setEntry(i, 1, getEntry(i, columnIndex));
            }
            return columnMatrix;
        }
    }

    /**
     * Get the transpose of the matrix (the matrix flipped on its diagonal)
     *
     * @return: the transpose of the matrix
     */
    Matrix getTranspose() {
        Matrix transp(columns, rows);
        // i and j start at 1 to stay consistent with matrix entry notation
        for (int i = 1; i <= transp.getRows(); i++) {
            for (int j = 1; j <= transp.getColumns(); j++) {
                transp.setEntry(i, j, getEntry(j, i));
            }
        }
        return transp;
    }

    /**
     * Check whether or not a matrix is a square matrix
     *
     * @return 1 (true) if the matrix is square, 0 (false) otherwise
     */
    bool isSquare() {
        return rows == columns;
    }

    /**
     * Check whether or not a matrix is a zero matrix
     *
     * @return 1 (true) if the matrix is a zero matrix, 0 (false) otherwise
     */
    bool isZero() {
        // check if each entry is 0
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < columns; j++) {
                if (matrix[i][j] != 0) {
                    return false;
                }
                else {
                    // pass
                }
            }
        }
        return true;
    }

    /**
     * Check whether or not a matrix is an identity matrix
     * (1's on the diagonal, 0's everywhere else)
     *
     * @return 1 (true) if the matrix is an identity matrix, 0 (false) otherwise
     */
    bool isIdentity() {
        if (isSquare()) {
            for (int i = 0; i < rows; i++) {
                for (int j = 0; j < columns; j++) {
                    // check for 1's on the diagonal
                    if (i == j) {
                        if (matrix[i][j] != 1) {
                            return false;
                        }
                        else {
                            // pass
                        }
                    }
                    // check for 0's everywhere else
                    else {
                        if (matrix[i][j] != 0) {
                            return false;
                        }
                        else {
                            // pass
                        }
                    }
                }
            }
            return true;
        }
        else {
            return false;
        }
    }

    /**
     * Turn the matrix into a zero matrix
     *
     * @postcond: the matrix will have zeros in all entries
     */
    void makeZero() {
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < columns; j++) {
                matrix[i][j] = 0;
            }
        }
    }

    /**
     * Turn the matrix into an identity matrix
     *
     * @precond: isSquare()
     * @postcond: if the matrix is a square matrix, the matrix will be turned into an identity matrix
     * (1's on the diagonal, 0's everywhere else)
     */
    void makeIdentity() {
        if (isSquare()) {
            for (int i = 0; i < rows; i++) {
                for (int j = 0; j < columns; j++) {
                    if (i == j) {
                        matrix[i][j] = 1;
                    }
                    else {
                        matrix[i][j] = 0;
                    }
                }
            }
        }
        else {
            throw 507;
        }
    }

    /**
     * Multiply all of the entries of the matrix by a scalar
     *
     * @postcond: all of the entries of the matrix are multiplied by s
     */
    void scalarMultiply(int s) {
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < columns; j++) {
                matrix[i][j] *= s;
            }
        }
    }

    /**
     * Compute the determinant of a matrix,
     * this method is static because it uses recursion,
     * and requires a matrix parameter in the recursive case
     *
     * @param m: the matrix to compute the determinant of
     * @precond: m.isSquare()
     * @return: the determinant of the matrix
     */
    static int determinant(Matrix m) {
        if (m.isSquare()) {
            // base case: zero matrix: return 0
            if (m.isZero()) {
                return 0;
            }
            // base case: identity matrix: return 1
            else if (m.isIdentity()) {
                return 1;
            }
            // base case: 1x1 matrix: return the lone entry
            else if (m.getRows() == 1) {
                return m.getEntry(1, 1);
            }
            // base case: 2x2 matrix
            else if (m.getRows() == 2) {
                return m.getEntry(1, 1) * m.getEntry(2, 2) - m.getEntry(1, 2) * m.getEntry(2, 1);
            }
            // recursive case: nxn matrix for n >= 3
            else {
                int det = 0;
                bool pastOmittedColumn = false;
                bool addNextDet = true;
                // i, j, k start at 1 to stay consistent with matrix entry notation
                for (int i = 1; i <= m.getColumns(); i++) {
                    // set up an n-1xn-1 matrix
                    Matrix subMatrix(m.getRows() - 1, m.getColumns() - 1);
                    for (int j = 1; j <= subMatrix.getRows(); j++) {
                        pastOmittedColumn = false;
                        for (int k = 1; k <= subMatrix.getColumns(); k++) {
                            // subMatrix will omit the first row and the ith column of the original matrix
                            if (i == k) {
                                pastOmittedColumn = true;
                            }
                            else {
                                // pass
                            }
                            if (!pastOmittedColumn) {
                                subMatrix.setEntry(j, k, m.getEntry(j+1, k));
                            }
                            else {
                                subMatrix.setEntry(j, k, m.getEntry(j+1, k+1));
                            }
                        }
                    }
                    // starting with adding,
                    // alternate between adding and subtractring the
                    // 1,i entry of original matrix times determinant of subMatrix
                    if (addNextDet) {
                        det += m.getEntry(1, i) * determinant(subMatrix);
                        addNextDet = false;
                    }
                    else {
                        det -= m.getEntry(1, i) * determinant(subMatrix);
                        addNextDet = true;
                    }
                    // deallocate the matrix before making a new one
                    subMatrix.deallocateMatrix();
                }
                return det;
            }
        }
        else {
            throw 508;
        }
    }

    /**
     * Check if two matrices have the same dimensions
     * (ex. if m1 is an mxn matrix, then m2 is also an mxn matrix)
     *
     * @param m1, m2: two matrices to be compared
     * @return: 1 (true) if the two matrices have the same number of rows and columns,
     * 0 (false) otherwise
     */
    static bool sameDimensions(Matrix m1, Matrix m2) {
        return m1.getRows() == m2.getRows() && m1.getColumns() == m2.getColumns();
    }

    /**
     * Check if two matrices are the same (i.e. their entries are the same)
     *
     * @param m1, m2: two matrices to be compared
     * @return: 1 (true) if the two matrices are the same,
     * 0 (false) otherwise
     */
    static bool sameMatrix(Matrix m1, Matrix m2) {
        if (sameDimensions(m1, m2)) {
            // check each entry of both matrices
            // i and j start at 1 to stay consistent with matrix entry notation
            for (int i = 1; i <= m1.getRows(); i++) {
                for (int j = 1; j <= m1.getColumns(); j++) {
                    if (m1.getEntry(i, j) != m2.getEntry(i, j)) {
                        return false;
                    }
                    else {
                        // pass
                    }
                }
            }
            return true;
        }
        else {
            return false;
        }
    }

    /**
     * Check if two matrices have the opposite dimensions
     * (ex. if m1 is an mxn matrix, then m2 is an nxm matrix)
     *
     * @param m1, m2: two matrices to be compared
     * @return: 1 (true) if the two matrices have the opposite dimensions,
     * 0 (false) otherwise
     */
    static bool oppositeDimensions(Matrix m1, Matrix m2) {
        return m1.getRows() == m2.getColumns() && m1.getColumns() == m2.getRows();
    }

    /**
     * Check if two matrices are transposes of each other
     *
     * @param m1, m2: two matrices to be compared
     * @return: 1 (true) if the two matrices are tranposes of each other,
     * 0 (false) otherwise
     */
    static bool transposes(Matrix m1, Matrix m2) {
        if (oppositeDimensions(m1, m2)) {
            // check each entry of both matrices
            // i and j start at 1 to stay consistent with matrix entry notation
            for (int i = 1; i <= m1.getRows(); i++) {
                for (int j = 1; j <= m1.getColumns(); j++) {
                    if (m1.getEntry(i, j) != m2.getEntry(j, i)) {
                        return false;
                    }
                    else {
                        // pass
                    }
                }
            }
            return true;
        }
        else {
            return false;
        }
    }

    /**
     * Check if a dot product can be computed from the given two matrices
     * (in this program, m1 has to be a 1xn matrix,
     * and m2 has to be an nx1 matrix)
     *
     * @param m1, m2: two matrices to check
     * @return: 1 (true) if a dot product can be computed from the two matrices,
     * 0 (false) otherwise
     */
    static bool canDot(Matrix m1, Matrix m2) {
        return oppositeDimensions(m1, m2) && m1.getRows() == 1 && m2.getColumns() == 1;
    }

    /**
     * Check if two matrices can be multiplied together
     * (in this program, if m1 is an mxn matrix and m2 is a pxq matrix,
     * then n must equal p)
     *
     * @param m1, m2: two matrices to check
     * @return: 1 (true) if the two matrices can be multiplied,
     * 0 (false) otherwise
     */
    static bool canMultiply(Matrix m1, Matrix m2) {
        return m1.getColumns() == m2.getRows();
    }

    /**
     * Add two matrices together and return their sum
     *
     * @param m1, m2: two matrices to be added
     * @precond: sameDimensions()
     * @return: a Matrix object that is the sum of m1 and m2
     */
    static Matrix sumOfMatrices(Matrix m1, Matrix m2) {
        if (sameDimensions(m1, m2)) {
            Matrix matrixSum(m1.getRows(), m1.getColumns());
            // i and j start at 1 to stay consistent with matrix entry notation
            for (int i = 1; i <= matrixSum.getRows(); i++) {
                for (int j = 1; j <= matrixSum.getColumns(); j++) {
                    matrixSum.setEntry(i, j, m1.getEntry(i, j) + m2.getEntry(i, j));
                }
            }
            return matrixSum;
        }
        else {
            throw 509;
        }
    }

    /**
     * Compute the dot product of two matrices
     *
     * @param m1, m2: two matrices to compute the dot product of
     * @precond: canDot()
     * @return: the dot product of the two matrices
     */
    static int dotProduct(Matrix m1, Matrix m2) {
        if (canDot(m1, m2)) {
            int dp = 0;
            // i starts at 1 to stay consistent with matrix entry notation
            for (int i = 1; i <= m1.getColumns(); i++) {
                dp += m1.getEntry(1, i) * m2.getEntry(i, 1);
            }
            return dp;
        }
        else {
            throw 510;
        }
    }

    /**
     * Compute the product of two matrices
     *
     * @param m1, m2: two matrices to compute the product of
     * @precond: canMultiply()
     * @return: a Matrix object that is the product of the two matrices,
     * (if m1 is an mxn matrix and m2 is a nxp matrix,
     * then this method returns an mxp matrix)
     */
    static Matrix productOfMatrices(Matrix m1, Matrix m2) {
        if (canMultiply(m1, m2)) {
            Matrix matrixProduct(m1.getRows(), m2.getColumns());
            // i and j start at 1 to stay consistent with matrix entry notation
            for (int i = 1; i <= matrixProduct.getRows(); i++) {
                Matrix rowMat = m1.getRowMatrix(i);
                for (int j = 1; j <= matrixProduct.getColumns(); j++) {
                    Matrix columnMat = m2.getColumnMatrix(j);
                    matrixProduct.setEntry(i, j, dotProduct(rowMat, columnMat));
                    // deallocate the column matrix before making a new one
                    columnMat.deallocateMatrix();
                }
                // deallocate the row matrix before making a new one
                rowMat.deallocateMatrix();
            }
            return matrixProduct;
        }
        else {
            throw 511;
        }
    }

    /**
     * Print the matrix to the console
     *
     * @postcond: the matrix is printed to the console
     */
    void printMatrix() {
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < columns; j++) {
                std::cout << matrix[i][j] << " ";
            }
            std::cout << "\n";
        }
    }

    /**
     * Print a matrix to the console
     *
     * @postcond: the specified matrix is printed to the console
     */
    static void staticPrintMatrix(Matrix m) {
        for (int i = 0; i < m.getRows(); i++) {
            for (int j = 0; j < m.getColumns(); j++) {
                std::cout << m.getMatrix()[i][j] << " ";
            }
            std::cout << "\n";
        }
    }

    /**
     * Create a deep clone of the matrix
     *
     * @return: a deep clone of the matrix
     */
    Matrix deepCloneMatrix() {
        Matrix mClone(rows, columns);
        // i and j start at 1 to stay consistent with matrix entry notation
        for (int i = 1; i <= rows; i++) {
            for (int j = 1; j <= columns; j++) {
                mClone.setEntry(i, j, getEntry(i, j));
            }
        }
        return mClone;
    }

    /**
     * Deallocate the memory allocated for the matrix when you are finished using it to prevent memory leaks,
     * DO NOT continue using the matrix after deallocating it
     *
     * @postcond: matrix will be deallocated and set to NULL
     */
    void deallocateMatrix() {
        // deallocate each inner array
        for (int i = 0; i < rows; i++) {
            delete[] matrix[i];
        }
        // deallocate the main array of pointers, and set matrix to NULL
        delete[] matrix;
        matrix = NULL;
    }
};

////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////// TEST SUITE //////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int main() {
    std::cout << "Test Suite for matrix.cpp\n";
    int uErrors = 0; // number of unintentional errors

    // create a couple of invalid Matrix objects
    std::cout << "Testing catching exceptions thrown by the Matrix constructor.\n";
    try {
        Matrix mInvalid1(0, 1);
        std::cout << "Error: Matrix() constructor did not throw an exception even though the number of rows was not positive.\n";
        uErrors += 1;
        mInvalid1.deallocateMatrix();
    }
    catch (int errNum) {
        // expected result
        std::cout << "Caught error number " << errNum << ".\n";
    }

    try {
        Matrix mInvalid2(1, 0);
        std::cout << "Error: Matrix() constructor did not throw an exception even though the number of columns was not positive.\n";
        uErrors += 1;
        mInvalid2.deallocateMatrix();
    }
    catch (int errNum) {
        // expected result
        std::cout << "Caught error number " << errNum << ".\n";
    }

    // create a valid Matrix object
    try {
        Matrix m1(2,3);

        // print the matrix
        std::cout << "Here is an example matrix.\n";
        m1.printMatrix();

        // test getRows()
        if (m1.getRows() != 2) {
            std::cout << "Error: getRows() returned " << m1.getRows() << " instead of 2.\n";
            uErrors += 1;
        }
        else {
            // expected result
        }

        // test getColumns()
        if (m1.getColumns() != 3) {
            std::cout << "Error: getColumns() returned " << m1.getColumns() << " instead of 3.\n";
            uErrors += 1;
        }
        else {
            // expected result
        }

        // test getMatrix()
        std::cout << "Here is the matrix being accessed and printed using getMatrix().\n";
        for (int i = 0; i < m1.getRows(); i++) {
            for (int j = 0; j < m1.getColumns(); j++) {
                std::cout << m1.getMatrix()[i][j] << " ";
            }
            std::cout << "\n";
        }

        // test staticPrintMatrix()
        std::cout << "Here is the matrix being printed using staticPrintMatrix().\n";
        Matrix::staticPrintMatrix(m1);

        // test IsSquare()
        if (m1.isSquare()) {
            std::cout << "Error: isSquare() returned true for non-square matrix.\n";
            uErrors += 1;
        }
        else {
            // expected result
        }

        // test IsZero()
        if (m1.isZero()) {
            // expected result
        }
        else {
            std::cout << "Error: isZero() returned false for zero matrix.\n";
            uErrors += 1;
        }

        // test IsIdentity()
        if (m1.isIdentity()) {
            std::cout << "Error: isIdentity() returned true for non-identity matrix.\n";
            uErrors += 1;
        }
        else {
            // expected result
        }

        // test sameDimensions()
        if (!Matrix::sameDimensions(m1, m1)) {
            std::cout << "Error: sameDimensions() returned false when comparing a matrix with itself.\n";
            uErrors += 1;
        }
        else {
            // expected result
        }

        // test sameMatrix()
        if (!Matrix::sameMatrix(m1, m1)) {
            std::cout << "Error: sameMatrix() returned false when comparing a matrix with itself.\n";
            uErrors += 1;
        }
        else {
            // expected result
        }

        // test oppositeDimensions()
        if (Matrix::oppositeDimensions(m1, m1)) {
            std::cout << "Error: oppositeDimensions() returned true when comparing a 3x2 matrix with itself.\n";
            uErrors += 1;
        }
        else {
            // expected result
        }

        // test transposes()
        if (Matrix::transposes(m1, m1)) {
            std::cout << "Error: transposes() returned true when comparing matrices that were not transposes.\n";
            uErrors += 1;
        }
        else {
            // expected result
        }

        // test setEntry()
        // change some entries
        try {
            m1.setEntry(1, 1, 1);
            m1.setEntry(1, 2, 2);
            m1.setEntry(2, 3, 3);
        }
        catch (int errNum) {
            std::cout << "Error: setEntry() threw an exception even though the specified entry was in bounds.\n";
            std::cout << "Error Number: " << errNum << ".\n";
            uErrors += 1;
        }
        // some invalid changes
        std::cout << "Testing catching exception thrown by setEntry().\n";
        try {
            m1.setEntry(0, 1, 7);
            std::cout << "Error: setEntry() did not throw an exception for an entry that was out of bounds.\n";
            uErrors += 1;
        }
        catch (int errNum) {
            // expected result
            std::cout << "Caught error number " << errNum << ".\n";
        }
        try {
            m1.setEntry(1, 0, 7);
            std::cout << "Error: setEntry() did not throw an exception for an entry that was out of bounds.\n";
            uErrors += 1;
        }
        catch (int errNum) {
            // expected result
        }
        try {
            m1.setEntry(3, 1, 7);
            std::cout << "Error: setEntry() did not throw an exception for an entry that was out of bounds.\n";
            uErrors += 1;
        }
        catch (int errNum) {
            // expected result
        }
        try {
            m1.setEntry(1, 4, 7);
            std::cout << "Error: setEntry() did not throw an exception for an entry that was out of bounds.\n";
            uErrors += 1;
        }
        catch (int errNum) {
            // expected result
        }
        
        // print the matrix
        std::cout << "After some changes...\n";
        m1.printMatrix();

        // test getEntry()
        std::cout << "The number at index 1,1 is " << m1.getEntry(1, 1) << ".\n";
        std::cout << "The number at index 2,2 is " << m1.getEntry(2, 2) << ".\n";
        std::cout << "The number at index 2,3 is " << m1.getEntry(2, 3) << ".\n";
        std::cout << "Testing catching exception thrown by getEntry().\n";
        try {
            m1.getEntry(0, 1);
            std::cout << "Error: getEntry() did not throw an exception for an entry that was out of bounds.\n";
            uErrors += 1;
        }
        catch (int errNum) {
            // expected result
            std::cout << "Caught error number " << errNum << ".\n";
        }
        try {
            m1.getEntry(1, 0);
            std::cout << "Error: getEntry() did not throw an exception for an entry that was out of bounds.\n";
            uErrors += 1;
        }
        catch (int errNum) {
            // expected result
        }
        try {
            m1.getEntry(3, 1);
            std::cout << "Error: getEntry() did not throw an exception for an entry that was out of bounds.\n";
            uErrors += 1;
        }
        catch (int errNum) {
            // expected result
        }
        try {
            m1.getEntry(1, 4);
            std::cout << "Error: getEntry() did not throw an exception for an entry that was out of bounds.\n";
            uErrors += 1;
        }
        catch (int errNum) {
            // expected result
        }

        // test IsZero()
        if (m1.isZero()) {
            std::cout << "Error: isZero() returned true for non-zero matrix.\n";
            uErrors += 1;
        }
        else {
            // expected result
        }

        // test scalarMultiply()
        m1.scalarMultiply(2);
        std::cout << "Scalar multiply the matrix by 2.\n";
        m1.printMatrix();

        // test makeZero()
        m1.makeZero();
        std::cout << "Turn matrix back into a zero matrix.\n";
        m1.printMatrix();

        // test IsZero()
        if (m1.isZero()) {
            // expected result
        }
        else {
            std::cout << "Error: isZero() returned false for zero matrix.\n";
            uErrors += 1;
        }

        // test makeIdentity()
        std::cout << "Testing catching exception thrown by makeIdentity().\n";
        try {
            m1.makeIdentity();
            std::cout << "Error: makeIdentity() did not throw an exception for a non-square matrix.\n";
            uErrors += 1;
        }
        catch (int errNum) {
            std::cout << "Caught error number " << errNum << ".\n";
            std::cout << "makeIdentity() will not work on a non-square matrix.\n";
        }

        // check if matrix is still a zero matrix
        std::cout << "Check if the matrix is still the zero matrix...\n";
        m1.printMatrix();

        // test IsZero()
        if (m1.isZero()) {
            // expected result
        }
        else {
            std::cout << "Error: isZero() returned false for zero matrix.\n";
            uErrors += 1;
        }

        // test determinant()
        std::cout << "Testing catching exception thrown by determinant().\n";
        try {
            int m1det = Matrix::determinant(m1);
            std::cout << "Error: determinant() did not throw an exception for a non-square matrix, but instead returned " << m1det << ".\n";
            uErrors += 1;
        }
        catch (int errNum) {
            std::cout << "Caught error number " << errNum << ".\n";
        }

        // deallocate the matrix when you are done using it
        // test deallocateMatrix()
        std::cout << "Deallocate the matrix at the end.\n";
        m1.deallocateMatrix();
        if (m1.getMatrix() != NULL) {
            std::cout << "Error: deallocateMatrix did not work because the matrix is not NULL.\n";
            uErrors += 1;
        }
        else {
            // expected result
        }
    }
    catch (int errNum) {
        std::cout << "Error: Matrix() constructor threw an exception even though both parameters were positive.\n";
        std::cout << "Error Number: " << errNum << ".\n";
        uErrors += 1;
    }

    // create a square matrix
    try {
        Matrix m2(4, 4);

        // print the matrix
        std::cout << "Here is a square matrix.\n";
        m2.printMatrix();

        // test getRows()
        if (m2.getRows() != 4) {
            std::cout << "Error: getRows() returned " << m2.getRows() << " instead of 4.\n";
            uErrors += 1;
        }
        else {
            // expected result
        }

        // test getColumns()
        if (m2.getColumns() != 4) {
            std::cout << "Error: getColumns() returned " << m2.getColumns() << " instead of 4.\n";
            uErrors += 1;
        }
        else {
            // expected result
        }

        // test getMatrix()
        std::cout << "Here is the matrix being accessed and printed using getMatrix().\n";
        for (int i = 0; i < m2.getRows(); i++) {
            for (int j = 0; j < m2.getColumns(); j++) {
                std::cout << m2.getMatrix()[i][j] << " ";
            }
            std::cout << "\n";
        }

        // test staticPrintMatrix()
        std::cout << "Here is the matrix being printed using staticPrintMatrix().\n";
        Matrix::staticPrintMatrix(m2);

        // test IsSquare()
        if (m2.isSquare()) {
            // expected result
        }
        else {
            std::cout << "Error: isSquare() returned false for square matrix.\n";
            uErrors += 1;
        }

        // test IsZero()
        if (m2.isZero()) {
            // expected result
        }
        else {
            std::cout << "Error: isZero() returned false for zero matrix.\n";
            uErrors += 1;
        }

        // test IsIdentity()
        if (m2.isIdentity()) {
            std::cout << "Error: isIdentity() returned true for non-identity matrix.\n";
            uErrors += 1;
        }
        else {
            // expected result
        }

        // test sameDimensions()
        if (!Matrix::sameDimensions(m2, m2)) {
            std::cout << "Error: sameDimensions() returned false when comparing a matrix with itself.\n";
            uErrors += 1;
        }
        else {
            // expected result
        }

        // test sameMatrix()
        if (!Matrix::sameMatrix(m2, m2)) {
            std::cout << "Error: sameMatrix() returned false when comparing a matrix with itself.\n";
            uErrors += 1;
        }
        else {
            // expected result
        }

        // test oppositeDimensions()
        if (!Matrix::oppositeDimensions(m2, m2)) {
            std::cout << "Error: oppositeDimensions() returned false when comparing a 4x4 matrix with itself.\n";
            uErrors += 1;
        }
        else {
            // expected result
        }

        // test transposes()
        if (!Matrix::transposes(m2, m2)) {
            std::cout << "Error: transposes() returned false when comparing matrices that were transposes.\n";
            uErrors += 1;
        }
        else {
            // expected result
        }

        // test setEntry()
        // change some entries
        try {
            m2.setEntry(1, 1, 4);
            m2.setEntry(3, 3, 2);
            m2.setEntry(4, 2, 7);
            m2.setEntry(1, 4, 1);
        }
        catch (int errNum) {
            std::cout << "Error: setEntry() threw an exception even though the specified entry was in bounds.\n";
            std::cout << "Error Number: " << errNum << ".\n";
            uErrors += 1;
        }
        // some invalid changes
        std::cout << "Testing catching exception thrown by setEntry().\n";
        try {
            m2.setEntry(0, 1, 8);
            std::cout << "Error: setEntry() did not throw an exception for an entry that was out of bounds.\n";
            uErrors += 1;
        }
        catch (int errNum) {
            // expected result
            std::cout << "Caught error number " << errNum << ".\n";
        }
        try {
            m2.setEntry(1, 0, 8);
            std::cout << "Error: setEntry() did not throw an exception for an entry that was out of bounds.\n";
            uErrors += 1;
        }
        catch (int errNum) {
            // expected result
        }
        try {
            m2.setEntry(5, 1, 8);
            std::cout << "Error: setEntry() did not throw an exception for an entry that was out of bounds.\n";
            uErrors += 1;
        }
        catch (int errNum) {
            // expected result
        }
        try {
            m2.setEntry(1, 5, 8);
            std::cout << "Error: setEntry() did not throw an exception for an entry that was out of bounds.\n";
            uErrors += 1;
        }
        catch (int errNum) {
            // expected result
        }
        
        // print the matrix
        std::cout << "After some changes...\n";
        m2.printMatrix();

        // test getEntry()
        std::cout << "The number at index 1,1 is " << m2.getEntry(1, 1) << ".\n";
        std::cout << "The number at index 4,2 is " << m2.getEntry(4, 2) << ".\n";
        std::cout << "The number at index 3,4 is " << m2.getEntry(3, 4) << ".\n";
        std::cout << "Testing catching exception thrown by getEntry().\n";
        try {
            m2.getEntry(0, 1);
            std::cout << "Error: getEntry() did not throw an exception for an entry that was out of bounds.\n";
            uErrors += 1;
        }
        catch (int errNum) {
            // expected result
            std::cout << "Caught error number " << errNum << ".\n";
        }
        try {
            m2.getEntry(1, 0);
            std::cout << "Error: getEntry() did not throw an exception for an entry that was out of bounds.\n";
            uErrors += 1;
        }
        catch (int errNum) {
            // expected result
        }
        try {
            m2.getEntry(5, 1);
            std::cout << "Error: getEntry() did not throw an exception for an entry that was out of bounds.\n";
            uErrors += 1;
        }
        catch (int errNum) {
            // expected result
        }
        try {
            m2.getEntry(1, 5);
            std::cout << "Error: getEntry() did not throw an exception for an entry that was out of bounds.\n";
            uErrors += 1;
        }
        catch (int errNum) {
            // expected result
        }

        // test IsZero()
        if (m2.isZero()) {
            std::cout << "Error: isZero() returned true for non-zero matrix.\n";
            uErrors += 1;
        }
        else {
            // expected result
        }

        // test scalarMultiply()
        m2.scalarMultiply(5);
        std::cout << "Scalar multiply the matrix by 5.\n";
        m2.printMatrix();

        // test transposes()
        if (Matrix::transposes(m2, m2)) {
            std::cout << "Error: transposes() returned true when comparing matrices that were not transposes.\n";
            uErrors += 1;
        }
        else {
            // expected result
        }

        // test makeIdentity()
        std::cout << "Turn the matrix into an identity matrix.\n";
        try {
            m2.makeIdentity();
            m2.printMatrix();
        }
        catch (int errNum) {
            std::cout << "Error: makeIdentity() threw an exception for a square matrix.\n";
            std::cout << "Error Number: " << errNum << ".\n";
            uErrors += 1;
        }
        
        // test isIdentity()
        if (m2.isIdentity()) {
            // expected result
        }
        else {
            std::cout << "Error: isIdentity() returned false for an identity matrix.\n";
            uErrors += 1;
        }

        // test IsZero()
        if (m2.isZero()) {
            std::cout << "Error: isZero() returned true for non-zero matrix.\n";
            uErrors += 1;
        }
        else {
            // expected result
        }

        // test determinant()
        try {
            int m2detIdentity = Matrix::determinant(m2);
            std::cout << "The determinant of this matrix is " << m2detIdentity << ".\n";
        }
        catch (int errNum) {
            std::cout << "Error: determinant() threw an exception for a square matrix.\n";
            std::cout << "Error Number: " << errNum << ".\n";
            uErrors += 1;
        }

        // test makeZero()
        m2.makeZero();
        std::cout << "Turn matrix back into a zero matrix.\n";
        m2.printMatrix();

        // test isIdentity()
        if (m2.isIdentity()) {
            std::cout << "Error: isIdentity() returned true for a zero matrix.\n";
            uErrors += 1;
        }
        else {
            // expected result
        }

        // test IsZero()
        if (m2.isZero()) {
            // expected result
        }
        else {
            std::cout << "Error: isZero() returned false for zero matrix.\n";
            uErrors += 1;
        }

        // test determinant()
        try {
            int m2detZero = Matrix::determinant(m2);
            std::cout << "The determinant of this matrix is " << m2detZero << ".\n";
        }
        catch (int errNum) {
            std::cout << "Error: determinant() threw an exception for a square matrix.\n";
            std::cout << "Error Number: " << errNum << ".\n";
            uErrors += 1;
        }

        // deallocate the matrix when you are done using it
        // test deallocateMatrix()
        std::cout << "Deallocate the matrix at the end.\n";
        m2.deallocateMatrix();
        // check if deallocateMatrix() worked
        if (m2.getMatrix() != NULL) {
            std::cout << "Error: deallocateMatrix did not work because the matrix is not NULL.\n";
            uErrors += 1;
        }
        else {
            // expected result
        }
    }
    catch (int errNum) {
        std::cout << "Error: Matrix() constructor threw an exception even though both parameters were positive.\n";
        std::cout << "Error Number: " << errNum << ".\n";
        uErrors += 1;
    }

    // test getRowMatrix() and getColumnMatrix
    Matrix m3(3,2);
    m3.setEntry(1,1,2);
    m3.setEntry(1,2,5);
    m3.setEntry(2,1,6);
    m3.setEntry(2,2,1);
    m3.setEntry(3,1,2);
    m3.setEntry(3,2,8);
    std::cout << "Here is an example matrix.\n";
    m3.printMatrix();
    try {
        Matrix m3row1 = m3.getRowMatrix(1);
        Matrix m3row2 = m3.getRowMatrix(2);
        Matrix m3row3 = m3.getRowMatrix(3);
        std::cout << "The first row of the matrix is\n";
        m3row1.printMatrix();
        std::cout << "The second row of the matrix is\n";
        m3row2.printMatrix();
        std::cout << "The third row of the matrix is\n";
        m3row3.printMatrix();
        // deallocate matrices
        m3row1.deallocateMatrix();
        m3row2.deallocateMatrix();
        m3row3.deallocateMatrix();
    }
    catch(int errNum) {
        std::cout << "Error: getRowMatrix() threw an exception for a valid row.\n";
        std::cout << "Error Number: " << errNum << ".\n";
        uErrors += 1;
    }
    try {
        Matrix m3col1 = m3.getColumnMatrix(1);
        Matrix m3col2 = m3.getColumnMatrix(2);
        std::cout << "The first column of the matrix is\n";
        m3col1.printMatrix();
        std::cout << "The second column of the matrix is\n";
        m3col2.printMatrix();
        // deallocate matrices
        m3col1.deallocateMatrix();
        m3col2.deallocateMatrix();
    }
    catch (int errNum) {
        std::cout << "Error: getColumnMatrix() threw an exception for a valid column.\n";
        std::cout << "Error Number: " << errNum << ".\n";
        uErrors += 1;
    }
    std::cout << "Testing catching exception thrown by getRowMatrix().\n";
    try {
        Matrix m3row0 = m3.getRowMatrix(0);
        std::cout << "Error: getRowMatrix() did not throw an exception for an invalid row.\n";
        uErrors += 1;
        m3row0.deallocateMatrix();
    }
    catch (int errNum) {
        // expected result
        std::cout << "Caught error number " << errNum << ".\n";
    }
    try {
        Matrix m3row4 = m3.getRowMatrix(4);
        std::cout << "Error: getRowMatrix() did not throw an exception for an invalid row.\n";
        uErrors += 1;
        m3row4.deallocateMatrix();
    }
    catch (int errNum) {
        // expected result
    }
    std::cout << "Testing catching exception thrown by getColumnMatrix().\n";
    try {
        Matrix m3col0 = m3.getColumnMatrix(0);
        std::cout << "Error: getColumnMatrix() did not throw an exception for an invalid column.\n";
        uErrors += 1;
        m3col0.deallocateMatrix();
    }
    catch (int errNum) {
        // expected result
        std::cout << "Caught error number " << errNum << ".\n";
    }
    try {
        Matrix m3col3 = m3.getColumnMatrix(3);
        std::cout << "Error: getColumnMatrix() did not throw an exception for an invalid column.\n";
        uErrors += 1;
        m3col3.deallocateMatrix();
    }
    catch (int errNum) {
        // expected result
    }

    // test deepCloneMatrix()
    try {
        Matrix m3Clone = m3.deepCloneMatrix();
        std::cout << "Here is a deep clone of the matrix.\n";
        m3Clone.printMatrix();

        // test sameDimensions()
        if (!Matrix::sameDimensions(m3, m3Clone)) {
            std::cout << "Error: sameDimensions() returned false for a matrix and its clone.\n";
            uErrors += 1;
        }

        // test sameMatrix()
        if (!Matrix::sameMatrix(m3, m3Clone)) {
            std::cout << "Error: sameMatrix() returned false for a matrix and its clone.\n";
            uErrors += 1;
        }

        // deallocate deep clone
        m3Clone.deallocateMatrix();
    }
    catch (int errNum) {
        std::cout << "Error: deepCloneMatrix() threw an exception.\n";
        std::cout << "Error Number: " << errNum << ".\n";
        uErrors += 1;
    }

    // test getTranspose()
    try {
        Matrix m3Transpose = m3.getTranspose();
        std::cout << "The matrix has this transpose.\n";
        m3Transpose.printMatrix();
        
        // test sameDimensions()
        if (Matrix::sameDimensions(m3, m3Transpose)) {
            std::cout << "Error: sameDimensions() returned true for transposes with different dimensions.\n";
            uErrors += 1;
        }
        else {
            // expected result
        }

        // test sameMatrix()
        if (Matrix::sameMatrix(m3, m3Transpose)) {
            std::cout << "Error: sameMatrix() returned true for transposes that are different.\n";
            uErrors += 1;
        }
        else {
            // expected result
        }

        // test oppositeDimensions()
        if (!Matrix::oppositeDimensions(m3, m3Transpose)) {
            std::cout << "Error: oppositeDiimensions() returned false for transposes.\n";
            uErrors += 1;
        }
        else {
            // expected result
        }

        // test transposes()
        if (!Matrix::transposes(m3, m3Transpose)) {
            std::cout << "Error: transposes() returned false for transposes.\n";
            uErrors += 1;
        }
        else {
            // expected result
        }

        // deallocate transpose
        m3Transpose.deallocateMatrix();
    }
    catch (int errNum) {
        std::cout << "Error: getTranspose() threw an exception.\n";
        std::cout << "Error Number: " << errNum << ".\n";
        uErrors += 1;
    }
    
    // deallocate the matrix when you are done using it
    std::cout << "Deallocate the matrix at the end.\n";
    m3.deallocateMatrix();

    // test determinant()
    Matrix m1x1(1, 1);
    Matrix m2x2(2, 2);
    Matrix m3x3(3, 3);
    Matrix m4x4(4, 4);
    
    m1x1.setEntry(1, 1, 5);
    
    m2x2.setEntry(1, 1, 3);
    m2x2.setEntry(1, 2, 5);
    m2x2.setEntry(2, 1, 2);
    m2x2.setEntry(2, 2, 4);

    m3x3.setEntry(1, 1, 7);
    m3x3.setEntry(1, 2, 5);
    m3x3.setEntry(1, 3, 3);
    m3x3.setEntry(2, 1, 4);
    m3x3.setEntry(2, 2, 3);
    m3x3.setEntry(2, 3, 2);
    m3x3.setEntry(3, 1, 6);
    m3x3.setEntry(3, 2, 9);
    m3x3.setEntry(3, 3, 9);

    m4x4.setEntry(1, 1, 1);
    m4x4.setEntry(1, 2, 2);
    m4x4.setEntry(1, 3, 4);
    m4x4.setEntry(1, 4, 8);
    m4x4.setEntry(2, 1, 9);
    m4x4.setEntry(2, 2, 7);
    m4x4.setEntry(2, 3, 5);
    m4x4.setEntry(2, 4, 3);
    m4x4.setEntry(3, 1, 6);
    m4x4.setEntry(3, 2, 4);
    m4x4.setEntry(3, 3, 3);
    m4x4.setEntry(3, 4, 2);
    m4x4.setEntry(4, 1, 1);
    m4x4.setEntry(4, 2, 6);
    m4x4.setEntry(4, 3, 9);
    m4x4.setEntry(4, 4, 9);

    std::cout << "Here are four matrices and their determinants.\n";
    int expectedDet;
    int returnedDet;
    
    m1x1.printMatrix();
    expectedDet = 5;
    try {
        returnedDet = Matrix::determinant(m1x1);
        std::cout << "Determinant: " << returnedDet << ".\n";
        if (returnedDet != expectedDet) {
            std::cout << "Error: determinant() returned the wrong determinant.\n";
            std::cout << "Correct determinant: " << expectedDet << ".\n";
            uErrors += 1;
        }
        else {
            // expected result
        }
    }
    catch (int errNum) {
        std::cout << "Error: determinant() threw an error for a square matrix.\n";
        std::cout << "Error Number: " << errNum << ".\n";
        uErrors += 1;
    }
    
    m2x2.printMatrix();
    expectedDet = 2;
    try {
        returnedDet = Matrix::determinant(m2x2);
        std::cout << "Determinant: " << returnedDet << ".\n";
        if (returnedDet != expectedDet) {
            std::cout << "Error: determinant() returned the wrong determinant.\n";
            std::cout << "Correct determinant: " << expectedDet << ".\n";
            uErrors += 1;
        }
        else {
            // expected result
        }
    }
    catch (int errNum) {
        std::cout << "Error: determinant() threw an error for a square matrix.\n";
        std::cout << "Error Number: " << errNum << ".\n";
        uErrors += 1;
    }
    
    m3x3.printMatrix();
    expectedDet = -3;
    try {
        returnedDet = Matrix::determinant(m3x3);
        std::cout << "Determinant: " << returnedDet << ".\n";
        if (returnedDet != expectedDet) {
            std::cout << "Error: determinant() returned the wrong determinant.\n";
            std::cout << "Correct determinant: " << expectedDet << ".\n";
            uErrors += 1;
        }
        else {
            // expected result
        }
    }
    catch (int errNum) {
        std::cout << "Error: determinant() threw an error for a square matrix.\n";
        std::cout << "Error Number: " << errNum << ".\n";
        uErrors += 1;
    }
    
    m4x4.printMatrix();
    expectedDet = 121;
    try {
        returnedDet = Matrix::determinant(m4x4);
        std::cout << "Determinant: " << returnedDet << ".\n";
        if (returnedDet != expectedDet) {
            std::cout << "Error: determinant() returned the wrong determinant.\n";
            std::cout << "Correct determinant: " << expectedDet << ".\n";
            uErrors += 1;
        }
        else {
            // expected result
        }
    }
    catch (int errNum) {
        std::cout << "Error: determinant() threw an error for a square matrix.\n";
        std::cout << "Error Number: " << errNum << ".\n";
        uErrors += 1;
    }

    // test canDot()
    if (Matrix::canDot(m4x4, m4x4)) {
        std::cout << "Error: canDot() returned true for two matrices that are not 1xn and nx1.\n";
        uErrors += 1;
    }
    else {
        // expected result
    }

    // deallocate matrices
    std::cout << "Deallocate the matrices at the end.\n";
    m1x1.deallocateMatrix();
    m2x2.deallocateMatrix();
    m3x3.deallocateMatrix();
    m4x4.deallocateMatrix();

    // test sumOfMatrices()
    Matrix m4(2, 3);
    Matrix m5(2, 3);
    Matrix mSumExpected(2, 3);

    m4.setEntry(1, 1, 1);
    m4.setEntry(1, 2, 2);
    m4.setEntry(1, 3, 3);
    m4.setEntry(2, 1, 4);
    m4.setEntry(2, 2, 5);
    m4.setEntry(2, 3, 6);

    m5.setEntry(1, 1, 7);
    m5.setEntry(1, 2, 8);
    m5.setEntry(1, 3, 9);
    m5.setEntry(2, 1, 1);
    m5.setEntry(2, 2, 0);
    m5.setEntry(2, 3, 1);

    mSumExpected.setEntry(1, 1, 8);
    mSumExpected.setEntry(1, 2, 10);
    mSumExpected.setEntry(1, 3, 12);
    mSumExpected.setEntry(2, 1, 5);
    mSumExpected.setEntry(2, 2, 5);
    mSumExpected.setEntry(2, 3, 7);

    std::cout << "Here are two matrices.\n";
    std::cout << "First matrix.\n";
    m4.printMatrix();
    std::cout << "Second matrix.\n";
    m5.printMatrix();
    std::cout << "Here is their sum.\n";
    try {
        Matrix mSum = Matrix::sumOfMatrices(m4, m5);
        mSum.printMatrix();
        if (!Matrix::sameMatrix(mSum, mSumExpected)) {
            std::cout << "Error: sumOfMatrices() returned the wrong sum.\n";
            std::cout << "Correct Sum:\n";
            mSumExpected.printMatrix();
            uErrors += 1;
        }
        else {
            // expected result
        }
        // deallocate sum
        mSum.deallocateMatrix();
    }
    catch (int errNum) {
        std::cout << "Error: sumOfMatrices() threw an exception for matrices with the same dimensions.\n";
        std::cout << "Error Number: " << errNum << ".\n";
        uErrors += 1;
    }

    // test canDot()
    if (Matrix::canDot(m4, m5)) {
        std::cout << "Error: canDot() returned true for two matrices that are not 1xn and nx1.\n";
        uErrors += 1;
    }
    else {
        // expected result
    }

    // test canMultiply()
    if (Matrix::canMultiply(m4, m5)) {
        std::cout << "Error: canMultiply() returned true for two matrices that cannot be multiplied.\n";
        uErrors += 1;
    }
    else {
        // expected result
    }

    // test dotProduct()
    std::cout << "Testing catching exception thrown by dotProduct().\n";
    try {
        int badDotProduct = Matrix::dotProduct(m4, m5);
        std::cout << "Error: dotProduct() did not throw an exception for matrices that cannot have their dot product calculated, but returned " << badDotProduct << ".\n";
        uErrors += 1;
    }
    catch (int errNum) {
        // expected result
        std::cout << "Caught error number " << errNum << ".\n";
    }

    // test productOfMatrices()
    std::cout << "Testing catching exception thrown by productOfMatrices().\n";
    try {
        Matrix badProduct = Matrix::productOfMatrices(m4, m5);
        std::cout << "Error: productOfMatrices() did not throw an exception for matrices that cannot be multiplied, but returned\n";
        badProduct.printMatrix();
        uErrors += 1;
        badProduct.deallocateMatrix();
    }
    catch (int errNum) {
        // expected result
        std::cout << "Caught error number " << errNum << ".\n";
    }

    // deallocate matrices
    std::cout << "Deallocate the matrices at the end.\n";
    m4.deallocateMatrix();
    m5.deallocateMatrix();
    mSumExpected.deallocateMatrix();

    // test canDot()
    Matrix mRow1(1, 3);
    Matrix mCol1(3, 1);
    Matrix mRow2(1, 2);
    Matrix mCol2(2, 1);

    mRow1.setEntry(1, 1, 6);
    mRow1.setEntry(1, 2, 4);
    mRow1.setEntry(1, 3, 2);

    mCol1.setEntry(1, 1, 2);
    mCol1.setEntry(2, 1, 5);
    mCol1.setEntry(3, 1, 0);

    mRow2.setEntry(1, 1, 7);
    mRow2.setEntry(1, 2, 3);

    mCol2.setEntry(1, 1, 8);
    mCol2.setEntry(2, 1, 1);

    if (Matrix::canDot(mRow1, mCol2)) {
        std::cout << "Error: canDot() returned true for two matrices that do not have opposite dimensions.\n";
        uErrors += 1;
    }
    else {
        // expected result
    }
    if (!Matrix::canDot(mRow1, mCol1)) {
        std::cout << "Error: canDot() returned false for two matrices that can have their dot product calculated.\n";
        uErrors += 1;
    }
    else {
        // expected result
    }
    if (!Matrix::canDot(mRow2, mCol2)) {
        std::cout << "Error: canDot() returned false for two matrices that can have their dot product calculated.\n";
        uErrors += 1;
    }
    else {
        // expected result
    }

    // test dotProduct()
    try {
        int badDotProduct = Matrix::dotProduct(mRow1, mCol2);
        std::cout << "Error: dotProduct() did not throw an exception for two matrices that do not have opposite dimensions, but returned " << badDotProduct << ".\n";
        uErrors += 1;
    }
    catch (int errNum) {
        // expected result
    }

    int expectedDotProduct;
    std::cout << "Here are two matrices.\n";
    std::cout << "First matrix.\n";
    mRow1.printMatrix();
    std::cout << "Second matrix.\n";
    mCol1.printMatrix();
    try {
        int dotProduct1 = Matrix::dotProduct(mRow1, mCol1);
        expectedDotProduct = 32;
        std::cout << "Their dot product is " << dotProduct1 << ".\n";
        if (dotProduct1 != expectedDotProduct) {
            std::cout << "Error: dotProduct() returned the wrong dot product.\n";
            std::cout << "Correct dot product: " << expectedDotProduct << ".\n";
            uErrors += 1;
        }
    }
    catch (int errNum) {
        std::cout << "Error: dotProduct() threw an exception for two matrices that can have their dot product calculated.\n";
        std::cout << "Error Number: " << errNum << ".\n";
        uErrors += 1;
    }

    std::cout << "Here are two matrices.\n";
    std::cout << "First matrix.\n";
    mRow2.printMatrix();
    std::cout << "Second matrix.\n";
    mCol2.printMatrix();
    try {
        int dotProduct2 = Matrix::dotProduct(mRow2, mCol2);
        expectedDotProduct = 59;
        std::cout << "Their dot product is " << dotProduct2 << ".\n";
        if (dotProduct2 != expectedDotProduct) {
            std::cout << "Error: dotProduct() returned the wrong dot product.\n";
            std::cout << "Correct dot product: " << expectedDotProduct << ".\n";
            uErrors += 1;
        }
    }
    catch (int errNum) {
        std::cout << "Error: dotProduct() threw an exception for two matrices that can have their dot product calculated.\n";
        std::cout << "Error Number: " << errNum << ".\n";
        uErrors += 1;
    }

    // test sumOfMatrices()
    std::cout << "Testing catching exception thrown by sumOfMatrices().\n";
    try {
        Matrix badSum = Matrix::sumOfMatrices(mRow1, mCol1);
        std::cout << "Error: sumOfMatrices() did not throw an exception for matrices of different dimensions, but returned\n";
        badSum.printMatrix();
        uErrors += 1;
        badSum.deallocateMatrix();
    }
    catch (int errNum) {
        // expected result
        std::cout << "Caught error number " << errNum << ".\n";
    }

    // test canMultiply()
    if (!Matrix::canMultiply(mRow1, mCol1)) {
        std::cout << "Error: canMultiply() returned false for two matrices that can be multiplied.\n";
        uErrors += 1;
    }
    else {
        // expected result
    }

    // deallocate matrices
    std::cout << "Deallocate the matrices at the end.\n";
    mRow1.deallocateMatrix();
    mCol1.deallocateMatrix();
    mRow2.deallocateMatrix();
    mCol2.deallocateMatrix();

    // test canMultiply()
    Matrix m6(2, 3);
    Matrix m7(3, 2);
    Matrix mProductExpected1(2, 2);
    Matrix mProductExpected2(3, 3);

    m6.setEntry(1, 1, 2);
    m6.setEntry(1, 2, 1);
    m6.setEntry(1, 3, 7);
    m6.setEntry(2, 1, 6);
    m6.setEntry(2, 2, 4);
    m6.setEntry(2, 3, 5);

    m7.setEntry(1, 1, 3);
    m7.setEntry(1, 2, 2);
    m7.setEntry(2, 1, 1);
    m7.setEntry(2, 2, 5);
    m7.setEntry(3, 1, 8);
    m7.setEntry(3, 2, 1);

    mProductExpected1.setEntry(1, 1, 63);
    mProductExpected1.setEntry(1, 2, 16);
    mProductExpected1.setEntry(2, 1, 62);
    mProductExpected1.setEntry(2, 2, 37);

    mProductExpected2.setEntry(1, 1, 18);
    mProductExpected2.setEntry(1, 2, 11);
    mProductExpected2.setEntry(1, 3, 31);
    mProductExpected2.setEntry(2, 1, 32);
    mProductExpected2.setEntry(2, 2, 21);
    mProductExpected2.setEntry(2, 3, 32);
    mProductExpected2.setEntry(3, 1, 22);
    mProductExpected2.setEntry(3, 2, 12);
    mProductExpected2.setEntry(3, 3, 61);

    if (!Matrix::canMultiply(m6, m7)) {
        std::cout << "Error: canMultiply() returned false for two matrices that can be multiplied.\n";
        uErrors += 1;
    }
    else {
        // expected result
    }
    if (!Matrix::canMultiply(m7, m6)) {
        std::cout << "Error: canMultiply() returned false for two matrices that can be multiplied.\n";
        uErrors += 1;
    }
    else {
        // expected result
    }

    // test productOfMatrices()
    std::cout << "Here are two matrices.\n";
    std::cout << "First matrix.\n";
    m6.printMatrix();
    std::cout << "Second matrix.\n";
    m7.printMatrix();
    try {
        Matrix mProduct1 = Matrix::productOfMatrices(m6, m7);
        std::cout << "The first matrix multiplied by the second matrix is\n";
        mProduct1.printMatrix();
        if (!Matrix::sameMatrix(mProduct1, mProductExpected1)) {
            std::cout << "Error: productOfMatrices() returned the wrong product.\n";
            std::cout << "Correct product: \n";
            mProductExpected1.printMatrix();
            uErrors += 1;
        }
        else {
            // expected result
        }
        // deallocate matrix
        mProduct1.deallocateMatrix();
    }
    catch (int errNum) {
        std::cout << "Error: productOfMatrices() threw an exception for two matrices that can be multiplied.\n";
        std::cout << "Error Number: " << errNum << ".\n";
        uErrors += 1;
    }
    try {
        Matrix mProduct2 = Matrix::productOfMatrices(m7, m6);
        std::cout << "The second matrix multiplied by the first matrix is\n";
        mProduct2.printMatrix();
        if (!Matrix::sameMatrix(mProduct2, mProductExpected2)) {
            std::cout << "Error: productOfMatrices() returned the wrong product.\n";
            std::cout << "Correct product: \n";
            mProductExpected2.printMatrix();
            uErrors += 1;
        }
        else {
            // expected result
        }
        // deallocate matrix
        mProduct2.deallocateMatrix();
    }
    catch (int errNum) {
        std::cout << "Error: productOfMatrices() threw an exception for two matrices that can be multiplied.\n";
        std::cout << "Error Number: " << errNum << ".\n";
        uErrors += 1;
    }

    // deallocate matrices
    std::cout << "Deallocate the matrices at the end.\n";
    m6.deallocateMatrix();
    m7.deallocateMatrix();
    mProductExpected1.deallocateMatrix();
    mProductExpected2.deallocateMatrix();

    std::cout << "Testing complete with " << uErrors << " unintentional errors.\n";
    return 0;
}
