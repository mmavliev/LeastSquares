#include <iostream>
#include <vector>
#include <iomanip>
#include <cmath>

using namespace std;

class Matrix {
public:
    int row, column;
    vector<vector<double>>matrix;
    Matrix() {
        row = 0, column = 0;
    }
    Matrix(int n, int m) {
        row = n;
        column = m;
        matrix.resize(n, vector<double>(m, 0));
    }

    friend istream& operator >> (istream& input, Matrix& cur_matrix) {
        input >> cur_matrix.row >> cur_matrix.column;
        cur_matrix.matrix.resize(cur_matrix.row, vector<double>(cur_matrix.column));
        for (int i = 0; i < cur_matrix.row; i++)
            for (int j = 0; j < cur_matrix.column; j++)
                input >> cur_matrix.matrix[i][j];
        return input;
    }

    friend ostream& operator << (ostream& output, Matrix& cur_matrix) {
        for (int i = 0; i < cur_matrix.row; i++) {
            for (int j = 0; j < cur_matrix.column; j++)
                output << cur_matrix.matrix[i][j] << " ";
            output << endl;
        }
        return output;
    }

    Matrix operator = (Matrix cur_matrix) {
        column = cur_matrix.column;
        row = cur_matrix.row;
        for (int i = 0; i < column; i++)
            for (int j = 0; j < row; j++)
                matrix[i][j] = cur_matrix.matrix[i][j];
        return *this;
    }

    friend Matrix operator + (Matrix first, Matrix second) {
        if (first.row != second.row || first.column != second.column)
            throw (string)("Error: the dimensional problem occurred\n");
        Matrix ans(first.row, first.column);
        for (int i = 0; i < first.row; i++)
            for (int j = 0; j < second.column; j++)
                ans.matrix[i][j] = first.matrix[i][j] + second.matrix[i][j];
        return ans;
    }

    friend Matrix operator - (Matrix first, Matrix second) {
        if (first.row != second.row || first.column != second.column)
            throw (string)("Error: the dimensional problem occurred\n");
        Matrix ans(first.row, first.column);
        for (int i = 0; i < first.row; i++)
            for (int j = 0; j < second.column; j++)
                ans.matrix[i][j] = first.matrix[i][j] - second.matrix[i][j];
        return ans;
    }

    friend Matrix operator * (Matrix first, Matrix second) {
        if (first.column != second.row)
            throw (string)("Error: the dimensional problem occurred\n");
        Matrix ans(first.row, second.column);
        for (int i = 0; i < first.row; i++)
            for (int j = 0; j < second.column; j++)
                for (int l = 0; l < first.matrix[i].size(); l++)
                    ans.matrix[i][j] += first.matrix[i][l] * second.matrix[l][j];
        return ans;
    }

    Matrix transpose() {
        Matrix ans(this->column, this->row);
        for (int i = 0; i < this->row; i++)
            for (int j = 0; j < this->column; j++)
                ans.matrix[j][i] = this->matrix[i][j];
        return ans;
    }
};

class SquareMatrix : public Matrix {
public:
    SquareMatrix(int x) {
        row = column = x;
        matrix.resize(x, vector<double>(x, 0));
    }
    SquareMatrix() {
        column = row = 0;
    }

    friend istream& operator >> (istream& input, SquareMatrix& cur_matrix) {
        input >> cur_matrix.row;
        cur_matrix.column = cur_matrix.row;
        cur_matrix.matrix.resize(cur_matrix.row, vector<double>(cur_matrix.row));
        for (int i = 0; i < cur_matrix.row; i++)
            for (int j = 0; j < cur_matrix.row; j++)
                input >> cur_matrix.matrix[i][j];
        return input;
    }

    friend ostream& operator << (ostream& output, SquareMatrix& cur_matrix) {
        for (int i = 0; i < cur_matrix.row; i++) {
            for (int j = 0; j < cur_matrix.row; j++)
                output << cur_matrix.matrix[i][j] << " ";
            output << endl;
        }
        return output;
    }

    friend SquareMatrix operator +(SquareMatrix first, SquareMatrix second) {
        Matrix* fir = &first, * sec = &second;
        Matrix ans = *fir + *sec;
        SquareMatrix* an = (SquareMatrix*)&ans;
        return *an;
    }

    friend SquareMatrix operator-(SquareMatrix first, SquareMatrix second) {
        Matrix* fir = &first, * sec = &second;
        Matrix ans = *fir - *sec;
        SquareMatrix* an = (SquareMatrix*)&ans;
        return *an;
    }

    friend SquareMatrix operator*(SquareMatrix first, SquareMatrix second) {
        Matrix* fir = &first, * sec = &second;
        Matrix ans = *fir * *sec;
        SquareMatrix* an = (SquareMatrix*)&ans;
        return *an;
    }

    SquareMatrix transpose() {
        Matrix fir = (*this);
        fir = (fir).transpose();
        SquareMatrix* ans = (SquareMatrix*)&fir;
        return *ans;
    }
    double GaussianDeterminant();
    void GaussianInverse();
private:
    int findMaxPivot(int i){
        int ind = i;
        double cur = matrix[i][i];
        for(int j = i; j < row; ++j){
            if(cur < abs(matrix[j][i])){
                cur = abs(matrix[j][i]);
                ind = j;
            }
        }
        return ind;
    }



};
Matrix unionMarix(SquareMatrix first, SquareMatrix second){
    Matrix resultOfUnion(first.row,first.column*2);
    for (int i=0;i<first.row;i++){
        for (int j=0;j<first.column;j++){
            resultOfUnion.matrix[i][j]=first.matrix[i][j];
            resultOfUnion.matrix[i][j+first.column]=second.matrix[i][j];
        }
    }
    return resultOfUnion;
}

class IdentityMatrix: public SquareMatrix{
public:
    IdentityMatrix(int dim):SquareMatrix(dim){
        for(int i = 0; i < dim; ++i){
            matrix[i][i] = 1;
        }
    }

};

class PermutationMatrix: public IdentityMatrix{
public:
    PermutationMatrix(Matrix& m, int f_row, int s_row): IdentityMatrix(m.row){
        swap(matrix[f_row], matrix[s_row]);
    }
};

class EliminationMatrix: public IdentityMatrix{
public:
    EliminationMatrix(Matrix& m, int mi_row, int sub_row): IdentityMatrix(m.row){
        matrix[mi_row][sub_row] = -(m.matrix[mi_row][sub_row]/m.matrix[sub_row][sub_row]);
    }
};

double SquareMatrix::GaussianDeterminant() {
    int step = 1;
    SquareMatrix result = *this;
    for(int i = 0; i < row; ++i){
        int ind = result.findMaxPivot(i);
        if(ind!=i){
            PermutationMatrix p(result, i, ind);
            result = p*result;
            cout << "step #"<<step++<<": permutation\n";
            cout <<fixed<<setprecision(2)<< result;
        }
        for(int k = i+1; k < row; ++k){
            if(result.matrix[k][i]==0.0){
                continue;
            } else {
                EliminationMatrix el(result, k, i);
                result = el*result;
                cout << "step #"<<step++<<": elimination\n";
                cout <<fixed<<setprecision(2)<< result;
            }
        }
    }
    double det = 1;
    for(int i = 0; i < row; ++i){
        det*= result.matrix[i][i];
    }
    return det;
}
void SquareMatrix::GaussianInverse() {
    //cout<<"step #0: Augmented Matrix\n";
    int step = 1;
    SquareMatrix result = *this;
    IdentityMatrix result1(result.row);
    SquareMatrix result2 = result1;
    Matrix newM(result.row,result.column*2);
    for (int i=0;i<result.row;i++){
        for (int j=0;j<result.column;j++){
            newM.matrix[i][j]=result.matrix[i][j];
            newM.matrix[i][j+result.column]=result2.matrix[i][j];
        }
    }
    //cout <<fixed<<setprecision(2)<<newM;
    //cout<<"Direct way:\n";

    for(int i = 0; i < row; ++i){
        int ind = result.findMaxPivot(i);
        if(ind!=i){
            PermutationMatrix p(result, i, ind);
            result = p*result;
            result2 = p*result2;
            Matrix newM(result.row,result.column*2);
            for (int i=0;i<result.row;i++){
                for (int j=0;j<result.column;j++){
                    newM.matrix[i][j]=result.matrix[i][j];
                    newM.matrix[i][j+result.column]=result2.matrix[i][j];
                }
            }

            //cout << "step #"<<step++<<": permutation\n";
            //cout <<fixed<<setprecision(2)<<newM;
        }
        for(int k = i+1; k < row; ++k){
            if(result.matrix[k][i]==0.0){
                continue;
            } else {
                EliminationMatrix el(result, k, i);
                result = el*result;
                result2 = el*result2;
                //cout << "step #"<<step++<<": elimination\n";
                Matrix newM(result.row,result.column*2);
                for (int i=0;i<result.row;i++) {
                    for (int j = 0; j < result.column; j++) {
                        newM.matrix[i][j] = result.matrix[i][j];
                        newM.matrix[i][j + result.column] = result2.matrix[i][j];
                    }
                }

                //cout <<fixed<<setprecision(2)<<newM;
            }
        }
    }
    //cout<<"Way back:\n";
    for(int i = row-1; i >0; i--){
        for(int k = i-1; k >=0; k--){
            if(result.matrix[k][i]==0.0){
                continue;
            } else {
                EliminationMatrix el(result, k, i);
                result = el*result;
                result2 = el*result2;
                //cout << "step #"<<step++<<": elimination\n";
                Matrix newM(result.row,result.column*2);
                for (int i=0;i<result.row;i++) {
                    for (int j = 0; j < result.column; j++) {
                        newM.matrix[i][j] = result.matrix[i][j];
                        newM.matrix[i][j + result.column] = result2.matrix[i][j];
                    }
                }
                //cout <<fixed<<setprecision(2)<< newM;
            }
        }

    }
    //cout<<"Diagonal normalization:\n";
    for (int i=0;i<result.row;i++){
        for (int j=0;j<result.column;j++){
            result2.matrix[i][j] = result2.matrix[i][j]/result.matrix[i][i];
        }
        result.matrix[i][i]=1.00;
    }
    Matrix newMwithNorm(result.row,result.column*2);
    for (int i=0;i<result.row;i++) {
        for (int j = 0; j < result.column; j++) {
            newMwithNorm.matrix[i][j] = result.matrix[i][j];
            newMwithNorm.matrix[i][j + result.column] = result2.matrix[i][j];
        }
    }
    //cout <<fixed<<setprecision(2)<< newMwithNorm;
    //cout<<"result:\n";
    //cout<<result2;
    *this = result2;
}

int main() {
    ios_base::sync_with_stdio(0);
    cin.tie(0);
    cout.tie(0);
    cout << fixed << setprecision(4);
    int m, n;
    cin >> m;
    vector<int> t(m);
    Matrix b(m, 1);
    for (int i = 0; i < m; i++) {
        cin >> t[i] >> b.matrix[i][0];
    }
    cin >> n;
    Matrix A(m, n+1);
    for (int i = 0; i < m; i++) {
        for (int j = 0; j <= n; j++) {
            A.matrix[i][j] = pow(t[i], j);
        }
    }
    cout << "A:" << endl;
    cout << A;
    cout << "A_T*A:" << endl;
    Matrix A1 = A.transpose()*A;
    cout << A1;
    cout << "(A_T*A)^-1:" << endl;
    SquareMatrix A2 = (*(SquareMatrix*)&A1);
    A2.GaussianInverse();
    cout << A2;
    cout << "A_T*b:" << endl;
    Matrix ATb = A.transpose()*b;
    cout << ATb;
    cout << "x~:" << endl;
    Matrix x = (Matrix)A2*ATb;
    cout << x;
	return 0;
}