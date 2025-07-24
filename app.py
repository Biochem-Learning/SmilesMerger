from flask import Flask, request, jsonify
from flask_cors import CORS
from process import combine_smile
import os

app = Flask(__name__)
CORS(app)

@app.route('/')
def index():
    return '<h1>Flask App is currently running</h1>'

@app.route('/combine_smile', methods=['POST'])
def combine_smile_route():
    data = request.json

    smile1 = data.get('smile1')
    smile2 = data.get('smile2')

    if not smile1 or not smile2:
        return jsonify({'error': 'Both SMILES must be provided'}), 400

    combined = combine_smile(smile1, smile2)
    return jsonify({'combined_smile': combined})

if __name__ == "__main__":
    app.run(debug=True, port=os.getenv("PORT", default=5000))
    