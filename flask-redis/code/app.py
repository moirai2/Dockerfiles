from flask import Flask
from redis import Redis

app = Flask(__name__, static_folder='static')
redis = Redis(host='redis', port=6379)

@app.route('/')
def hello():
    count = redis.incr('hits')
    return 'Hello World! I have been seen {} times. <img src="static/flask-logo.png">\n'.format(count)

if __name__ == "__main__":
    app.run(host="0.0.0.0", debug=True)